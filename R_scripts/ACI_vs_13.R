# Preamble #####################################################################

## For Jonathan

# file_path <- "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/"

# Packages #####################################################################

## Load the following packages
# library(readr)
# library(DiceKriging)

# Importing the inputs  ########################################################

## The normalised input settings from the PPE
# inputs_x_norm_T_matrix <- 
#   as.matrix(
#     read.csv(paste0(file_path, 
#                     "Design/ACURE-UKESM1-PPE-par-norm-Design.csv"
#                     )
#              )
#   )
# inputs_x_norm_matrix <- t(inputs_x_norm_T_matrix)

## There are p variables in all, and n runs
# p <- ncol(inputs_x_norm_T_matrix)
# n <- nrow(inputs_x_norm_T_matrix)

# Importing the outputs ########################################################

## Hd_oct is the observable variable
response_Hd_oct_dataframe <- 
  read.table(
    paste0(file_path, 
           "Optimal-Constraint-Data/Data_for_emulation/oct/CDNC_PPE_NH_SH_difference_oct_revised_match_MODIS_CF.dat"
           ), 
  col.names = "Hd_oct"
)

## ACI is the unobserved variable
# response_ACI_dataframe <- 
#   read.table(
#     paste0(file_path, 
#            "./Optimal-Constraint-Data/data_for_emulation/jul/ACI_PPE_global_jul.dat"
#     ), 
#     col.names = "ACI"
#   )

# lm work ######################################################################

## In case linear models are fine, and there's no need to fit GPs, we check 
# adjusted R^2 values

## For Hd_oct
R2adj_Hd_oct <- summary(
  lm(Hd_oct ~., data.frame(cbind(response_Hd_oct_dataframe, inputs_x_norm_T_matrix)))
)$adj.r.squared

## For ACI
# R2adj_ACI <- summary(
#   lm(ACI ~., data.frame(cbind("ACI" = response_ACI_dataframe, inputs_x_norm_T_matrix)))
# )$adj.r.squared

## We export both objects to /objects, with the import code below
write_lines(R2adj_Hd_oct,
            file="objects/R2adj_Hd_oct.txt")
# R2adj_Hd_oct <-
#   as.numeric(readLines("objects/R2adj_Hd_oct.txt"))
# write_lines(R2adj_ACI,
#             file="objects/R2adj_ACI.txt")
# R2adj_ACI <-
#   as.numeric(readLines("objects/R2adj_ACI.txt"))

## Both have low adjusted R^2 values, so we will fit GPs

# GP set-up ####################################################################

## We specify some variables, namely N (the size of the sample from the input
## parameter space) and q (the number of betas that will be estimated) 
# N <- 500
# q <- 1 + p

## We load the required packages for sampling
# library(lhs)

## We create the sample
# xstar_matrix <- t(randomLHS(N, p))
## Export this to /objects, with the import code below
# write_lines(as.vector(xstar_matrix),
#             file="objects/xstar_matrix.txt")
# xstar_matrix <-
#   matrix(as.numeric(readLines("objects/xstar_matrix.txt")),
#          nrow = p)

## Our estimated posterior derivative requires H, which doesn't depend on the 
# response
# H_matrix <- cbind(c(rep(1, n)), t(unname(inputs_x_norm_matrix)))

## We define a function called gp_d_dxi_normalised_function which returns the
## normalised GP estimates of the partial derivatives. It takes as its single
## argument a 221-long vector or 221 x 1 dataframe of responses

gp_d_dxi_normalised_function <- function(response, loglik_index) {

  gp <- km(~.,
           design = inputs_x_norm_T_matrix,
           response = response,
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
  )

  ## Using the posterior estimates of the range parameters, we make the matrix
  ## A, its inverse, and t(xstar)^T
  lhat_vector <- gp@covariance@range.val
  lhat_matrix <- matrix(lhat_vector, nrow = p)
  lhatinv_diag_matrix <- diag(1/lhat_vector)

  A_matrix <- matrix(NA, nrow = n, ncol = n)
  for (r in 1:n) {
    for (c in 1:n) {
      if (is.na(A_matrix[r,c])) {
        A_matrix[r,c] <- A_matrix[c,r] <-
          exp(-0.5 *
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                (lhatinv_diag_matrix^2) %*%
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
          )
      }
    }
  }
  Ainv_matrix <- solve(A_matrix)

  txstar_T_matrix <- matrix(NA, nrow = N, ncol = n)
  for (r in 1:N) {
    for (c in 1:n) {
      txstar_T_matrix[r,c] <- exp(-0.5 *
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                                    (lhatinv_diag_matrix^2) %*%
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
      )
    }
    cat("txstarT: ", r, "/", N, ",", sep = "")
  }

  ## Secondly, we use the trend estimates to estimate the 37 posterior partial
  ## derivatives at each of the N points, then normalise them
  betahat_matrix <- matrix(gp@trend.coef, nrow = q)

  pd_T_dataframe <- data.frame(matrix(NA, nrow = p, ncol = N))

  for (r in 1:p) {
    pd_T_dataframe[r,] <-
      matrix(betahat_matrix[r+1,], nrow = N) +
      (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (2*lhat_vector[r]^2)) *
         txstar_T_matrix) %*%
      Ainv_matrix %*% (as.matrix(response, ncol = n) - H_matrix %*% betahat_matrix)
    cat("pd_t:", r, "/", p, ",", sep = "")
  }

  pd_dataframe <- t(pd_T_dataframe)

  pdnormed_dataframe <- pd_dataframe / sqrt(rowSums(pd_dataframe^2))

  pdnormed_T_matrix <- t(data.matrix(pdnormed_dataframe))

  ## just when repeating runs
  # assign(paste0("loglik_long142lat88_", loglik_index), logLik.km(gp),
  #        envir = .GlobalEnv)

  return(pdnormed_T_matrix)
}

## First we'll fit the GP for ACI and calculate the normalised partial
## derivatives
gppdnormalised_t_ACI_matrix <- gp_d_dxi_normalised_function(response_ACI_dataframe)

## Next we'll fit a GPs for Hd_oct

gppdnormalised_t_Hd_oct_matrix <- gp_d_dxi_normalised_function(response_Hd_oct_dataframe)

# And now the alignment measure

AM_with_corr <- mean(abs(colSums(gppdnormalised_t_ACI_matrix * gppdnormalised_t_Hd_oct_matrix)))

gp_d_dxi_normalised_function <- function(response, loglik_index) {
  
  gp <- km(~.,
           design = inputs_x_norm_T_matrix,
           response = response,
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
  )
  
  ## Using the posterior estimates of the range parameters, we make the matrix
  ## A, its inverse, and t(xstar)^T
  lhat_vector <- gp@covariance@range.val
  lhat_matrix <- matrix(lhat_vector, nrow = p)
  lhatinv_diag_matrix <- diag(1/lhat_vector)
  
  A_matrix <- matrix(NA, nrow = n, ncol = n)
  for (r in 1:n) {
    for (c in 1:n) {
      if (is.na(A_matrix[r,c])) {
        A_matrix[r,c] <- A_matrix[c,r] <-
          exp(-0.5 *
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                (lhatinv_diag_matrix^2) %*%
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
          )
      }
    }
  }
  Ainv_matrix <- solve(A_matrix)
  
  txstar_T_matrix <- matrix(NA, nrow = N, ncol = n)
  for (r in 1:N) {
    for (c in 1:n) {
      txstar_T_matrix[r,c] <- exp(-0.5 *
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                                    (lhatinv_diag_matrix^2) %*%
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
      )
    }
    cat("txstarT: ", r, "/", N, ",", sep = "")
  }
  
  ## Secondly, we use the trend estimates to estimate the 37 posterior partial
  ## derivatives at each of the N points, then normalise them
  betahat_matrix <- matrix(gp@trend.coef, nrow = q)
  
  pd_T_dataframe <- data.frame(matrix(NA, nrow = p, ncol = N))
  
  for (r in 1:p) {
    pd_T_dataframe[r,] <-
      matrix(betahat_matrix[r+1,], nrow = N) +
      (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (lhat_vector[r]^2)) *
         txstar_T_matrix) %*%
      Ainv_matrix %*% (as.matrix(response, ncol = n) - H_matrix %*% betahat_matrix)
    cat("pd_t:", r, "/", p, ",", sep = "")
  }
  
  pd_dataframe <- t(pd_T_dataframe)
  
  pdnormed_dataframe <- pd_dataframe / sqrt(rowSums(pd_dataframe^2))
  
  pdnormed_T_matrix <- t(data.matrix(pdnormed_dataframe))
  
  ## just when repeating runs
  # assign(paste0("loglik_long142lat88_", loglik_index), logLik.km(gp),
  #        envir = .GlobalEnv)
  
  return(pdnormed_T_matrix)
}

## First we'll fit the GP for ACI and calculate the normalised partial
## derivatives
gppdnormalised_t_ACI_matrix <- gp_d_dxi_normalised_function(response_ACI_dataframe)

## Next we'll fit a GPs for Hd_oct

gppdnormalised_t_Hd_oct_matrix <- gp_d_dxi_normalised_function(response_Hd_oct_dataframe)

# And now the alignment measure

AM_wo_corr <- mean(abs(colSums(gppdnormalised_t_ACI_matrix * gppdnormalised_t_Hd_oct_matrix)))
