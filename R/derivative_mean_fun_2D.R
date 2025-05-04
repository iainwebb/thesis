# ------------
# a 2D example
# ------------

# SETUP
# number of observations ------------------
n <- 5
# specify the inputs ----------------------
p <- 2
x_1_vector <- seq(-.5,2*pi,len=n+2)[c(1,2,3,4,6)]
x_2_vector <- sample((seq(-.5,2*pi,len=n+2)[c(1,2,3,4,6)]), 5)
x_matrix <- matrix(c(x_1_vector, x_2_vector), nrow=p, ncol=n, byrow = T)
x_matrix_T <- t(x_matrix)
# real-life function ----------------------
real_life_process <- function(x){
  sin(x)
}
# observations ----------------------------
f_obs_vector <- x_1_vector + 3*x_2_vector
f_obs_matrix <- matrix(f_obs_vector, nrow=1)
# data ------------------------------------
f_data_dataframe <- data.frame(x=t(x_matrix),
                               f_obs=t(f_obs_matrix))

# PRIOR MEAN FUNCTION
# h function form ----------------------
## for zero mean, must be 1 + 0*x
q <- 3
h_function_form_part_1 <- function(x) {
  1 + 0*x
}
h_function_form_part_2 <- function(x_1, x_2) {
  1*x_1 + 1*x_2
}
# and corresponding formula
## for zero mean,doesn't matter if ~0 or ~1
h_formula <- ~.

# GP
# fit the GP ------------------------------
library(DiceKriging)

f_GP <- km(formula = h_formula, 
           design = data.frame(x = t(x_matrix)), 
           response = data.frame(f_obs_vector), 
           covtype = "gauss"
           )
f_GP
# extract estimates
betas_hat_matrix <- matrix(f_GP@trend.coef,
                    nrow = q)
sigma_sq_hat <- f_GP@covariance@sd2
l_hat_vector <- f_GP@covariance@range.val
l_hat_matrix <- matrix(f_GP@covariance@range.val,
       nrow = p)
l_hat_diag_matrix <- matrix(0, nrow = p, ncol = p)
l_hat_diag_matrix[1,1] <- l_hat_vector[1]
l_hat_diag_matrix[2,2] <- l_hat_vector[2]

# PLOT
# define x_star  ---------------------------
x_star_lower <- 0
x_star_upper <- 1
x.1_star_vector <- seq(from=x_star_lower, to=x_star_upper, length.out = 3)
x.2_star_vector <- seq(from=x_star_lower, to=x_star_upper, length.out = 3)

library(tidyr)
x_star_dataframe <- x_star_dataframe %>% expand(x.1, x.2)

# predict at x_star values using GP --------

x_star_predictions_list <- predict(f_GP,
                                     newdata = x_star_dataframe,
                                     type="SK"
                                     )
# make the plot ----------------------------
plot(x_star_vector, 
     x_star_predictions_list$mean, 
     type = "l", 
     ylim = c(-2,2),
     # ylim = c(min(x_star_predictions_list$lower95),max(x_star_predictions_list$upper95)),
     xlab = "x", 
     ylab = "y"
     )
abline(h=0)
lines(x_star_vector, x_star_predictions_list$trend, col="violet", lty=2, lwd=3)
# lines(x_star_vector, x_star_predictions_list$lower95, col=color$UK, lty=2)
# lines(x_star_vector, x_star_predictions_list$upper95, col=color$UK, lty=2)
points(x_vector, f_obs_vector, col="red", pch=19)

# RECREATING POSTERIOR MEAN FUNCTION
# defining covariance matrix making function
sq_exp_cov_function <- function(matrix_1, matrix_2, l = l_hat_diag_matrix){
  cov_matrix <- matrix(rep(0, ncol(matrix_1) * ncol(matrix_2)), nrow = ncol(matrix_1))
  for (a in 1:nrow(cov_matrix)) {
    for (b in 1:ncol(cov_matrix)) {
      cov_matrix[a,b] <- exp(-matrix((matrix_1[,a] - matrix_2[,b]), nrow = 1) %*% l_hat_diag_matrix %*% matrix((matrix_1[,a] - matrix_2[,b]), ncol = 1))
    }
  }
  return(cov_matrix)
}
# create A
A <- sq_exp_cov_function(x_matrix, x_matrix) + diag(0.001, ncol(x_matrix))
# make h(x_star)^T
h_x_star_T_matrix <- matrix(
  c(rep(1, nrow(x_star_dataframe)),
    x_star_dataframe[,1],
    x_star_dataframe[,2]
    ), 
  nrow = nrow(x_star_dataframe),
  ncol = q, 
  byrow = F)
# make t(x_star)^T
t_x_star_T_matrix <- sq_exp_cov_function(t(as.matrix(x_star_dataframe)), x_matrix)
# make H
H <- matrix(
  c(rep(1, ncol(x_matrix)),
    x_matrix[1,],
    x_matrix[2,]
  ), 
  nrow = ncol(x_matrix),
  ncol = q, 
  byrow = F)
# make predictions at x_star
x_star_predictions_manual_matrix <- h_x_star_T_matrix %*% betas_hat_matrix +  t_x_star_T_matrix %*% solve(A) %*% (f_data_dataframe$f_obs - H %*% betas_hat_matrix)
# check they are the same (difference = 0)
abs(x_star_predictions_list$mean - x_star_predictions_manual_matrix) > 0.00000000000001

# make the plot ----------------------------
plot(x_star_vector, 
     x_star_predictions_manual_matrix, 
     type = "l",
     ylim = c(-2,2),
     # ylim = c(min(x_star_predictions_list$lower95),max(x_star_predictions_list$upper95)),
     xlab = "x", 
     ylab = "y"
)
# check with predictions made with DiceKriging
# display every m points
m <- 50
points(x_star_vector[seq(1,length(x_star_vector),m)],
       x_star_predictions_list$mean[seq(1,length(x_star_vector),m)],
       pch = 4)
abline(h=0)
lines(x_star_vector, x_star_predictions_list$trend, col="violet", lty=2, lwd=3)
# lines(x_star_vector, x_star_predictions_list$lower95, col=color$UK, lty=2)
# lines(x_star_vector, x_star_predictions_list$upper95, col=color$UK, lty=2)
points(x_vector, f_obs_vector, col="red", pch=19)

# DERIVATIVE OF POSTERIOR MEAN
# we have betas_hat_matrix already
betas_hat_matrix
# have h(x_star)^T already
h_x_star_T_matrix
# get h ready to be differentiated ----------------------
h_expression_part_1 <- expression(1 + 0*x)
h_expression_part_2 <- expression(1*x)
# the derivative wrt x_i is 
d_dx_k_h_part_1 <- D(h_expression_part_1, 'x')
d_dx_k_h_part_2 <- D(h_expression_part_2, 'x')
# make d_dx_k_h_x_star
d_dx_k_h_x_star_matrix <- matrix(c(rep(d_dx_k_h_part_1, length(x_star_vector)), rep(d_dx_k_h_part_2, length(x_star_vector))), nrow = length(x_star_vector), ncol = q)

d_dx_i_h_x_star_T_matrix <- cbind(matrix(0, nrow = p), diag(1, p))

# have A_inv already
A_inv <- solve(A)
# have t(x_star)^T already
t_x_star_T_matrix
# make d/dx_i [t(x_star)^T]
N <- nrow(x_star_dataframe)

# make d_dx_i_m_star_m_star_matrix
i <- 2
k <- 8
j <- 4
x_star_dataframe

partial_derivatives_dataframe <- data.frame(d_dx.1 = rep(NA, length(x.1_star_vector) * length(x.2_star_vector)))
for (i in 2:p) {
  partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe, 
                                         data.frame(rep(NA, length(x.1_star_vector) * length(x.2_star_vector))))
}
colnames(partial_derivatives_dataframe) <- c("d_dx.1", "d_dx.2")

for (i in 1:p) {
  for (k in 1:N) {
    d_dx_i_t_x_star_T_matrix <- matrix(NA, ncol = n)
    for (j in 1:n) {
      d_dx_i_t_x_star_T_matrix[,j] <- -2 / l_hat_matrix[i,]^2 * (x_star_dataframe[[k,i]] - x_matrix[i,j]) * t_x_star_T_matrix[k,j]
    }
    partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_dataframe[k,i] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix %*% A_inv %*% (f_data_dataframe$f_obs - H %*% betas_hat_matrix)))
  }
}

partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))

partial_derivatives_normalised_response_surface_1_dataframe <- partial_derivatives_normalised_dataframe
partial_derivatives_normalised_response_surface_2_dataframe <- partial_derivatives_normalised_dataframe

abs(sum(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe)) / nrow(x_star_dataframe)

d_dx_i_h_x_star_T_matrix[i,] %*% betas_hat_matrix - 2/l_hat_matrix[i,]

d_dx_i_t_x_star_T_matrix <-
  matrix(
    NA, nrow = length(x_star_vector), ncol = n
  )
for (j in 1:n) {
  for (k in 1:N) {
    d_dx_i_t_x_star_T_matrix[k,j] <- 
      -1/l_hat^2* (x_star_vector[k] - x_vector[j]) * sq_exp_cov_function(x_star_vector[k], x_vector[j])
  }
}
# make d_dx_i_m_star_m_star_matrix
x_star_derivative_predictions_manual_matrix <- 
  d_dx_k_h_x_star_matrix %*% betas_hat_matrix + 
  d_dx_i_t_x_star_T_matrix %*% solve(A) %*% (f_data_dataframe$f_obs - H %*% betas_hat_matrix)
# make the plot ----------------------------
plot(x_star_vector, 
     x_star_derivative_predictions_manual_matrix, 
     type = "l",
     ylim = c(-2,2),
     # ylim = c(min(x_star_predictions_list$lower95),max(x_star_predictions_list$upper95)),
     xlab = "x", 
     ylab = "y"
)
# check with predictions made with DiceKriging
# display every m points
m <- 50
points(x_star_vector[seq(1,length(x_star_vector),m)],
       x_star_predictions_list$mean[seq(1,length(x_star_vector),m)],
       pch = 4)
abline(h=0)
lines(x_star_vector, x_star_predictions_list$trend, col="violet", lty=2, lwd=3)
# lines(x_star_vector, x_star_predictions_list$lower95, col=color$UK, lty=2)
# lines(x_star_vector, x_star_predictions_list$upper95, col=color$UK, lty=2)
points(x_vector, f_obs_vector, col="red", pch=19)
