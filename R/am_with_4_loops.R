rm(list = ls())
calcs <- 5

library(ncdf4)
library(DiceKriging)
library(lhs)
  
JDO_jul <- nc_open(paste("C:/Users/smp22ijw/Downloads/ACURE_P3_AOD_Total_jul.nc"))
JDO_jul <- ncvar_get(JDO_jul, names(JDO_jul$var))
AODs_Total_jul_level2 <- JDO_jul[,,2,]

# heatmap(t(AODs_Total_jul_level2[,,125]), Rowv = NA, Colv = NA)

# SETUP
# specify the inputs ----------------------
x_norm_T_matrix <- read.csv("C:/Users/smp22ijw/Downloads/ACURE-UKESM1-PPE-par-norm-Design.csv")
x_norm_matrix <- t(x_norm_T_matrix)
# x_norm_T_matrix <- data.frame(matrix(runif(221*37), nrow = 221))
colnames(x_norm_T_matrix) <- paste0(rep("x.",37), 1:37)
# x_norm_T_matrix <- x_norm_T_matrix[,1:15]
p <- ncol(x_norm_T_matrix)
# x_norm_T_matrix <- read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")

# number of observations ------------------
n <- nrow(x_norm_T_matrix)

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
x_star_T_dataframe <- data.frame(maximinLHS(N, p))
colnames(x_star_T_dataframe) <- paste0(rep("x.",p), 1:p)
x_star_matrix <- t(as.matrix(x_star_T_dataframe))
# make H_matrix
H_matrix <- matrix(
  c(rep(1, ncol(x_norm_matrix)),
    x_norm_matrix[1,], x_norm_matrix[2,], x_norm_matrix[3,], x_norm_matrix[4,], x_norm_matrix[5,], x_norm_matrix[6,], x_norm_matrix[7,], 
    x_norm_matrix[8,], x_norm_matrix[9,], x_norm_matrix[10,], x_norm_matrix[11,], x_norm_matrix[12,], x_norm_matrix[13,],
    x_norm_matrix[14,], x_norm_matrix[15,], x_norm_matrix[16,], x_norm_matrix[17,], x_norm_matrix[18,], x_norm_matrix[19,], 
    x_norm_matrix[20,], x_norm_matrix[21,], x_norm_matrix[22,], x_norm_matrix[23,], x_norm_matrix[24,], x_norm_matrix[25,], 
    x_norm_matrix[26,], x_norm_matrix[27,], x_norm_matrix[28,], x_norm_matrix[29,], x_norm_matrix[30,], x_norm_matrix[31,], 
    x_norm_matrix[32,], x_norm_matrix[33,], x_norm_matrix[34,], x_norm_matrix[35,], x_norm_matrix[36,], x_norm_matrix[37,]
  ), 
  nrow = ncol(x_norm_matrix),
  ncol = q, 
  byrow = F)

gridboxes <- matrix(c(100, 85, 100, 85), nrow = 2, byrow = T)
# gridboxes <- matrix(c(100, 85, 20, 120), nrow = 2, byrow = T)



for (c in 1:calcs) {
  
  start_time <- Sys.time()
  
  for (g in 1:2) {
    # SETUP
    # observations ----------------------------
    AODs_Total_jul_level2_gb <- AODs_Total_jul_level2[gridboxes[g,1], gridboxes[g,2],]
    
    # GP
    # fit the GP ------------------------------
    f_GP <- km(~., 
               design = x_norm_T_matrix, 
               response = AODs_Total_jul_level2_gb, 
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
                                       newdata = x_star_T_dataframe,
                                       type="SK"
    )
    # create A_matrix
    A_matrix <- sq_exp_cov_function(x_norm_matrix, x_norm_matrix) + diag(0.001, nrow(x_norm_T_matrix))
    # make A_inv_matrix
    A_inv_matrix <- solve(A_matrix)
    
    # partial derivatives
    partial_derivatives_dataframe <- data.frame(rep(NA, N))
    for (i in 2:p) {
      partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe, 
                                             data.frame(rep(NA, N)))
    }
    colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx.",37), 1:37)
    
    # make t(x_star)^T
    t_x_star_T_matrix <- sq_exp_cov_function(x_star_matrix, x_norm_matrix)
    
    for (i in 1:p) {
      for (k in 1:N) {
        d_dx_i_t_x_star_T_matrix <- matrix(NA, ncol = n)
        for (j in 1:n) {
          d_dx_i_t_x_star_T_matrix[,j] <- -2 / l_hat_matrix[i,]^2 * (x_star_T_dataframe[[k,i]] - x_norm_matrix[i,j]) * t_x_star_T_matrix[k,j]
        }
        partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_T_dataframe[k,i] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix %*% A_inv_matrix %*% (AODs_Total_jul_level2_gb - H_matrix %*% betas_hat_matrix)))
      }
    }
    
    partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
    assign(paste0("partial_derivatives_normalised_response_surface_", g, "_dataframe"), partial_derivatives_normalised_dataframe)
    
  }
  
AM <- abs(sum(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe)) / N
  
assign(paste0("AM_with_", N, "_attempt_", c), AM)
  
end_time <- Sys.time()
  
assign(paste0("time_with_", N, "_attempt_", c), end_time - start_time)
}

for (c in 1:calcs) {
  print(eval(parse(text = paste0("AM_with_", N, "_attempt_", c))))
}

for (c in 1:calcs) {
  print(eval(parse(text = paste0("time_with_", N, "_attempt_", c))))
}