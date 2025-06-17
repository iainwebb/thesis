# HIDE -------
# i <- 36
# k <- 9
# for (i in 1:p) {
#   for (k in 1:N) {
#     d_dx_i_t_x_star_T_matrix <- matrix(-2 / l_hat_matrix[i,]^2 * (rep(x_star_matrix[i,k],n) - x_norm_matrix[i,]) * t_x_star_T_matrix[k,], ncol = n)
#     partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix %*% A_inv_matrix %*% (AODs_Total_jul_level2_gb - H_matrix %*% betas_hat_matrix)))
#   }
# }
# 
# d_dx_i_t_x_star_T_matrix <- matrix(-2 / l_hat_matrix[i,]^2 * (rep(x_star_matrix[i,k],n) - x_norm_matrix[i,]) * t_x_star_T_matrix[k,], ncol = n)
# 
# 
# partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe
# 
# dim(abs(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe))

# START ------
setwd("U:/ManWin/My Documents/thesis")
rm(list = ls())
calcs <- 1

library(ncdf4)
library(DiceKriging)
library(lhs)

JDO_jul <- nc_open(paste("data/ACURE_P3_AOD_Total_jul.nc"))
JDO_jul <- ncvar_get(JDO_jul, names(JDO_jul$var))
AODs_Total_jul_level2 <- JDO_jul[,,2,]

# heatmap(t(AODs_Total_jul_level2[,,125]), Rowv = NA, Colv = NA)

# SETUP
# specify the inputs ----------------------
# x_norm_T_matrix <- read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")
x_norm_T_matrix <- read.csv("data/ACURE-UKESM1-PPE-par-norm-Design.csv")
#                                                                       [,1:4]
x_norm_matrix <- t(x_norm_T_matrix)
p <- ncol(x_norm_T_matrix)
colnames(x_norm_T_matrix) <- paste0(rep("x.",p), 1:p)

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
N <- 100000
x_star_matrix <- t(randomLHS(N, p))
# x_star_matrix <- t(read.csv("data/x_norm_matrix.csv", header = FALSE, sep = ""))[,sample(1:200000, 10000)]
N <- ncol(x_star_matrix)
x_star_T_dataframe <- data.frame(t(x_star_matrix))
colnames(x_star_T_dataframe) <- paste0(rep("x.",p), 1:p)
# make H_matrix
H_matrix <- cbind(c(rep(1, n)), t(unname(x_norm_matrix)))

gridboxes <- matrix(c(20, 120, 20, 120), nrow = 2, byrow = T)
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
    assign(paste0("betas_with_", N, "_gb_", g, "_attempt_", c), f_GP@trend.coef)
    assign(paste0("l_hat_with_", N, "_gb_", g, "_attempt_", c), f_GP@covariance@range.val)
    
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
    colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx.",p), 1:p)
    
    # make t(x_star)^T
    t_x_star_T_matrix <- sq_exp_cov_function(x_star_matrix, x_norm_matrix)
    
    # i <- 36
    # k <- 9
    for (i in 1:p) {
      d_dx_i_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
      for (k in 1:N) {
        partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (AODs_Total_jul_level2_gb - H_matrix %*% betas_hat_matrix)))
      }
    }
    
    partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
    assign(paste0("partial_derivatives_normalised_response_surface_", g, "_dataframe"), partial_derivatives_normalised_dataframe)
    
  }
  
  AM <- sum(abs(rowSums(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe))) / N
  
  # assign(paste0("AM_with_", N, "_attempt_", c), AM)
  write(AM, 
        file=paste0("R/land_with_itself_", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_AM_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(AM, 
        file=paste0("R/land_with_itself_", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("betas_with_", N, "_gb_", 1, "_attempt_", c)))), 
        file=paste0("R/land_with_itself_", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("betas_with_", N, "_gb_", 2, "_attempt_", c)))), 
        file=paste0("R/land_with_itself_", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 1, "_attempt_", c)))), 
        file=paste0("R/land_with_itself_", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 2, "_attempt_", c)))), 
        file=paste0("R/land_with_itself_", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  
  end_time <- Sys.time()
  
  # assign(paste0("time_with_", N, "_attempt_", c), end_time - start_time)
  
}

# for (c in 1:calcs) {
#   print(eval(parse(text = paste0("AM_with_", N, "_attempt_", c))))
#   print(eval(parse(text = paste0("betas_with_", N, "_gb_", 1, "_attempt_", c))))
#   print(eval(parse(text = paste0("betas_with_", N, "_gb_", 2, "_attempt_", c))))
#   print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 1, "_attempt_", c))))
#   print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 2, "_attempt_", c))))
# }
# 
# for (c in 1:calcs) {
#   print(eval(parse(text = paste0("time_with_", N, "_attempt_", c))))
# }
# 
# print(c(p, N))
