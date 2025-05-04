# ------------
# a 1D example
# ------------

# SETUP
# number of observations ------------------
n <- 5
# specify the inputs ----------------------
x_vector <- seq(-.5,2*pi,len=n+2)[c(1,2,3,4,6)]
x_matrix <- matrix(x_vector, ncol=1)
# real-life function ----------------------
real_life_process <- function(x){
  sin(x)
}
# observations ----------------------------
f_obs_vector <- real_life_process(x_matrix)
f_obs_matrix <- as.matrix(f_obs_vector, ncol=1)
# data ------------------------------------
f_data_dataframe <- data.frame(x=x_matrix,
                               f_obs=f_obs_matrix)

# PRIOR MEAN FUNCTION
# h function form ----------------------
## for zero mean, must be 1 + 0*x
q <- 2
h_function_form_part_1 <- function(x) {
  1 + 0*x
}
h_function_form_part_2 <- function(x) {
  1*x # 
}
# and corresponding formula
## for zero mean,doesn't matter if ~0 or ~1
h_formula <- ~.

# GP
# fit the GP ------------------------------
library(DiceKriging)
f_GP <- km(formula = h_formula, 
           design = data.frame(x = x_vector), 
           response = data.frame(f_obs_vector), 
           covtype = "gauss"
           )
# extract estimates
beta_0_hat <- f_GP@trend.coef[1]
beta_1_hat <- f_GP@trend.coef[2]
betas_hat_matrix <- matrix(c(beta_0_hat, beta_1_hat),
                    nrow = q)
sigma_sq_hat <- f_GP@covariance@sd2
l_hat <- f_GP@covariance@range.val

# PLOT
# define x_star  ---------------------------
x_star_lower <- 0 - 3
x_star_upper <- 2*pi + 3
x_star_vector <- seq(from=x_star_lower, to=x_star_upper, by=0.005)
x_star_dataframe <- data.frame(x = x_star_vector)
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
sq_exp_cov_function <- function(matrix_1, matrix_2, l = l_hat){
  cov_matrix <- matrix(rep(0, length(matrix_1) * length(matrix_2)), nrow = length(matrix_1))
  for (i in 1:nrow(cov_matrix)) {
    for (j in 1:ncol(cov_matrix)) {
      cov_matrix[i,j] <- exp(-(sqrt((matrix_1[i] - matrix_2[j])^2)/l)^2/2)
    }
  }
  return(cov_matrix)
}
# create A
A <- sq_exp_cov_function(x_vector, x_vector) + diag(0.001, length(x_vector))
# make h(x_star)^T
h_x_star_T_matrix <- matrix(
  c(h_function_form_part_1(x_star_vector),
    h_function_form_part_2(x_star_vector)), 
  nrow = length(x_star_vector),
  ncol = q, 
  byrow = F)
# make t(x_star)^T
t_x_star_T_matrix <- sq_exp_cov_function(x_star_vector, x_vector)
# make H
H <- matrix(c(h_function_form_part_1(x_vector),
              h_function_form_part_2(x_vector)), 
              nrow = length(x_vector),
              ncol = q)
# make predictions at x_star
x_star_predictions_manual_matrix <- h_x_star_T_matrix %*% betas_hat_matrix +  t_x_star_T_matrix %*% solve(A) %*% (f_data_dataframe$f_obs - H %*% betas_hat_matrix)
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
# get h ready to be diffentiated ----------------------
h_expression_part_1 <- expression(1 + 0*x)
h_expression_part_2 <- expression(1*x)
# the derivative wrt x_i is 
d_dx_k_h_part_1 <- D(h_expression_part_1, 'x')
d_dx_k_h_part_2 <- D(h_expression_part_2, 'x')
# make d_dx_k_h_x_star
d_dx_k_h_x_star_matrix <- matrix(c(rep(d_dx_k_h_part_1, length(x_star_vector)), rep(d_dx_k_h_part_2, length(x_star_vector))), nrow = length(x_star_vector), ncol = q)
# have A_inv already
A_inv <- solve(A)
# have t(x_star)^T already
t_x_star_T_matrix
# make d/dx_i [t(x_star)^T]
N <- length(x_star_vector)
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
