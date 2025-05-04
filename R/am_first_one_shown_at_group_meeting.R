rm(list = ls())

for (c in 1:5) {
  
start_time <- Sys.time()
  
  
# library(ncdf4)
# library(DiceKriging)
# library(lhs)

JDO_jul <- nc_open(paste("C:/Users/smp22ijw/Downloads/ACURE_P3_AOD_Total_jul.nc"))
JDO_jul <- ncvar_get(JDO_jul, names(JDO_jul$var))
AODs_Total_jul_level2 <- JDO_jul[,,2,]

# heatmap(t(AODs_Total_jul_level2[,,125]), Rowv = NA, Colv = NA)

# SETUP
# specify the inputs ----------------------
X_norm_T <- read.csv("C:/Users/smp22ijw/Downloads/ACURE-UKESM1-PPE-par-norm-Design.csv")
X_norm <- t(X_norm_T)
# X_norm_T <- data.frame(matrix(runif(221*37), nrow = 221))
colnames(X_norm_T) <- paste0(rep("x.",37), 1:37)
# X_norm_T <- X_norm_T[,1:15]
p <- ncol(X_norm_T)
# X_norm_T <- read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")

# number of observations ------------------
n <- nrow(X_norm_T)

# PRIOR MEAN FUNCTION
q <- 1 + p
# and corresponding formula
h_formula <- ~.

# DERIVATIVE WORK
# defining covariance matrix making function
sq_exp_cov_function <- function(matrix_1, matrix_2, l = l_hat_diag_matrix){
  cov_matrix <- matrix(rep(NA, ncol(matrix_1) * ncol(matrix_2)), nrow = ncol(matrix_1))
  for (a in 1:nrow(cov_matrix)) {
    for (b in 1:ncol(cov_matrix)) {
      cov_matrix[a,b] <- exp(-matrix((matrix_1[,a] - matrix_2[,b]), nrow = 1) %*% l_hat_diag_matrix %*% matrix((matrix_1[,a] - matrix_2[,b]), ncol = 1))
    }
  }
  return(cov_matrix)
}
# define x_star  ---------------------------
x_star_lower <- 0
x_star_upper <- 1
N <- 100
x_star_dataframe <- data.frame(maximinLHS(N, p))
colnames(x_star_dataframe) <- paste0(rep("x.",p), 1:p)
# make H
H <- matrix(
  c(rep(1, ncol(X_norm)),
    X_norm[1,], X_norm[2,], X_norm[3,], X_norm[4,], X_norm[5,], X_norm[6,], X_norm[7,], 
    X_norm[8,], X_norm[9,], X_norm[10,], X_norm[11,], X_norm[12,], X_norm[13,],
    X_norm[14,], X_norm[15,], X_norm[16,], X_norm[17,], X_norm[18,], X_norm[19,], 
    X_norm[20,], X_norm[21,], X_norm[22,], X_norm[23,], X_norm[24,], X_norm[25,], 
    X_norm[26,], X_norm[27,], X_norm[28,], X_norm[29,], X_norm[30,], X_norm[31,], 
    X_norm[32,], X_norm[33,], X_norm[34,], X_norm[35,], X_norm[36,], X_norm[37,]
  ), 
  nrow = ncol(X_norm),
  ncol = q, 
  byrow = F)

gridboxes <- matrix(c(100, 85, 100, 85), nrow = 2, byrow = T)
# gridboxes <- matrix(c(100, 85, 20, 120), nrow = 2, byrow = T)


start_time <- Sys.time()

for (g in 1:2) {
  # SETUP
  # observations ----------------------------
  AODs_Total_jul_level2_gb1 <- AODs_Total_jul_level2[gridboxes[g,1], gridboxes[g,2],]
  
  # GP
  # fit the GP ------------------------------
  f_GP <- km(~., 
             design = X_norm_T, 
             response = AODs_Total_jul_level2_gb1, 
             covtype="gauss", optim.method="BFGS", control=list(maxit=500))
  # extract estimates
  betas_hat_matrix <- matrix(f_GP@trend.coef,
                             nrow = q)
  sigma_sq_hat <- f_GP@covariance@sd2
  l_hat_vector <- f_GP@covariance@range.val
  l_hat_matrix <- matrix(f_GP@covariance@range.val,
                         nrow = p)
  l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
  
  # DERIVATIVE WORK
  # predict at x_star values using GP --------
  x_star_predictions_list <- predict(f_GP,
                                     newdata = x_star_dataframe,
                                     type="SK"
  )
  # create A
  A <- sq_exp_cov_function(X_norm, X_norm) + diag(0.001, nrow(X_norm_T))
  # make A_inv
  A_inv <- solve(A)
  
  # partial derivatives
  partial_derivatives_dataframe <- data.frame(rep(NA, N))
  for (i in 2:p) {
    partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe, 
                                           data.frame(rep(NA, N)))
  }
  colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx.",37), 1:37)
  
  # make t(x_star)^T
  t_x_star_T_matrix <- sq_exp_cov_function(t(as.matrix(x_star_dataframe)), X_norm)
  
  for (i in 1:p) {
    for (k in 1:N) {
      d_dx_i_t_x_star_T_matrix <- matrix(NA, ncol = n)
      for (j in 1:n) {
        d_dx_i_t_x_star_T_matrix[,j] <- -2 / l_hat_matrix[i,]^2 * (x_star_dataframe[[k,i]] - X_norm[i,j]) * t_x_star_T_matrix[k,j]
      }
      partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_dataframe[k,i] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix %*% A_inv %*% (AODs_Total_jul_level2_gb1 - H %*% betas_hat_matrix)))
    }
  }
  
  partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
  assign(paste0("partial_derivatives_normalised_response_surface_", g, "_dataframe"), partial_derivatives_normalised_dataframe)
  
}

AM <- abs(sum(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe)) / N

assign(paste0("AM_with_", N, "_attempt_", c), AM)

end_time <- Sys.time()

assign(paste0("time_with_", N, "_attempt_", c), cat(end_time - start_time))

}

AM_with_100_attempt_5
