#------------------------
# a 1 dimensional example
#------------------------
# dimensional of the inputs
dim_inputs <- 1    
# number of the inputs
num_obs <- 5
# uniform samples of design
input <- matrix(seq(0.1, 0.9, length.out=num_obs), num_obs)

# Following codes use maximin Latin Hypercube Design, which is typically better than uniform
# library(lhs)
# input <- maximinLHS(n=num_obs, k=dim_inputs)  ##maximin lhd sample

####
# outputs from 2sin(3x)
test_function_1D <- function(x) {
  2*sin(3*x - 0.5)
}
output <- matrix(test_function_1D(input), num_obs)
plot(output ~ input, xlim = c(0,1), pch=19, col = "blue")

# use constant mean basis, with no constraint on optimization
m1<- rgasp(design = input, response = output, lower_bound=FALSE,
           kernel_type = 'pow_exp', alpha = 2)

# the following use constraints on optimization
# m1<- rgasp(design = input, response = output, lower_bound=TRUE)

# the following use a single start on optimization
# m1<- rgasp(design = input, response = output, lower_bound=FALSE)

# number of points to be predicted 
num_testing_input <- 1000 
# generate points to be predicted
testing_input <- matrix(sort(runif(num_testing_input*dim_inputs)), num_testing_input,dim_inputs)
# Perform prediction
m1.predict<-predict(m1, testing_input, outasS3 = FALSE)
# Predictive mean
m1.predict@mean

## Comparing matern52 (default) to Gauss
m1_mat52<- rgasp(design = input, response = output, lower_bound=FALSE)
m1_mat52.predict <- predict(m1_mat52, testing_input, outasS3 = FALSE)
sum(abs(m1_mat52.predict@mean - m1.predict@mean))

# The following tests how good the prediction is 
testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-test_function_1D(testing_input[i,])
}

# compute the MSE, average coverage and average length
# out of sample MSE
MSE_emulator <- sum((m1.predict@mean-testing_output)^2)/(num_testing_input)

# proportion covered by 95% posterior predictive credible interval
prop_emulator <- length(which((m1.predict@lower95<=testing_output)
                              &(m1.predict@upper95>=testing_output)))/num_testing_input

plot(testing_input, test_function_1D(testing_input),
     main="1D RobustGaSP example",
     xlab="input",
     ylab="output",
     type="l",
     col="black")
lines(testing_input,m1.predict@mean, col="red")
lines(testing_input,m1.predict@lower95, col="red", lty = "dashed")
lines(testing_input,m1.predict@upper95, col="red", lty = "dashed")
points(input, output, pch=19, col="blue")
legend("topleft",
       c("actual","prediction","observations"),
       fill=c("black","red","blue")
)
plot(testing_input, test_function_1D(testing_input),
     main="1D RobustGaSP example",
     xlab="input",
     ylab="output",
     type="l",
     col="black")
lines(testing_input,m1_mat52.predict@mean, col="red")
lines(testing_input,m1_mat52.predict@lower95, col="red", lty = "dashed")
lines(testing_input,m1_mat52.predict@upper95, col="red", lty = "dashed")
points(input, output, pch=19, col="blue")
legend("topleft",
       c("actual","prediction","observations"),
       fill=c("black","red","blue")
)

# average length of posterior predictive credible interval
length_emulator <- sum(m1.predict@upper95-m1.predict@lower95)/num_testing_input

# output of prediction
MSE_emulator
prop_emulator
length_emulator  
# normalized RMSE
sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))

# checking against DiceKriging ----
# library("DiceKriging")
gp <- km(~., 
         input,
         output, 
         covtype="gauss", optim.method="BFGS", control=list(maxit=500)
)
# Perform prediction
gp.predict<-predict(gp, testing_input, "SK", checkNames = F)
# Predictive mean
gp.predict$mean

sum(abs(m1.predict@mean - gp.predict$mean)) / num_testing_input

# The following tests how good the prediction is 
testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-test_function_1D(testing_input[i,])
}

# compute the MSE, average coverage and average length
# out of sample MSE
MSE_emulator_km <- sum((gp.predict$mean-testing_output)^2)/(num_testing_input)

# proportion covered by 95% posterior predictive credible interval
prop_emulator_km <- length(which((gp.predict$lower95<=testing_output)
                              &(gp.predict$upper95>=testing_output)))/num_testing_input
prop_emulator_km

plot(testing_input, test_function_1D(testing_input),
     main="1D DiceKriging example",
     xlab="input",
     ylab="output",
     type="l",
     col="black")
lines(testing_input,gp.predict$mean, col="red")
lines(testing_input,gp.predict$lower95, col="red", lty = "dashed")
lines(testing_input,gp.predict$upper95, col="red", lty = "dashed")
points(input, output, pch=19, col="blue")
legend("topleft",
       c("actual","prediction","observations"),
       fill=c("black","red","blue")
)

## AM calculation ----

# Our estimated posterior derivative requires H
H_matrix <- cbind(c(rep(1, num_obs)), input)

# Comparing estimates
m1@

m1@sigma2_hat
gp@covariance@sd2

# We extract the length parameter estimates

# We make the matrix A
A_matrix <- matrix(NA, nrow = num_obs, ncol = num_obs)
r <- 4; c <- 2
for (r in 1:n) {
  for (c in 1:n) {
    if (is.na(A_matrix[r,c])) {
      A_matrix[r,c] <- A_matrix[c,r] <- 
        exp(-0.5*
              matrix(input[r,] - input[c,], nrow = 1) %*% ## MIGHT NEED AMENDING SINCE NOW WORKING WITH THE TRANSPOSE
              (lhatinv_diag_matrix^2) %*%
              matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
        )
    }
  }
}


## Using the posterior estimates of the range parameters, we make the matrix 
## A, its inverse, and t(xstar)^T
lhat_vector <- gp@covariance@range.val
lhat_matrix <- matrix(lhat_vector)
lhatinv_diag_matrix <- diag(x = 1/lhat_vector, nrow = p)
sqrt(gp@covariance@sd2)
sqrt(2)






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

betahat_matrix <- matrix(gp@trend.coef, nrow = q)

plot(response_test_dataframe$y ~ as.vector(inputs_x_norm_matrix), xlim = c(0,1), ylim = c(-2,2))

predictions <- predict(gp, as.vector(t(xstar_matrix)), "SK")
points(predictions$mean ~ as.vector(xstar_matrix), pch=19, cex=0.5)

my_predictions <- matrix(c(rep(1, N), xstar_matrix), nrow = N, byrow = F) %*% betahat_matrix +
  txstar_T_matrix %*% Ainv_matrix %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
points(my_predictions ~ as.vector(xstar_matrix[1,]), pch=19, cex=0.5, col = "red")

sum(abs(predictions$mean - as.vector(my_predictions))) / (10^(-8)*N)
plot((predictions$mean - as.vector(my_predictions)) ~ as.vector(xstar_matrix))
points(rep(0,5) ~ as.vector(inputs_x_norm_matrix), pch = 19, col = "red")

## Secondly, we use the trend estimates to estimate the 37 posterior partial
## derivatives at each of the N points, then normalise them

pd_T_dataframe <- data.frame(matrix(NA, nrow = p, ncol = N))

for (r in 1:p) {
  pd_T_dataframe[r,] <-
    matrix(betahat_matrix[r+1,], nrow = N) + 
    (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (lhat_vector[r]^2)) *
       txstar_T_matrix) %*% 
    Ainv_matrix %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
  cat("pd_t:", r, "/", p, ",", sep = "")
}


plot(response_test_dataframe$y ~ as.vector(inputs_x_norm_matrix), xlim = c(0,1), ylim = c(-20,20))

points((6*cos(as.vector(3*xstar_matrix-0.5))) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "blue")
points(as.numeric(as.vector(pd_T_dataframe[1,])) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "red")

as.numeric(as.vector(pd_T_dataframe[1,])) / (6*cos(as.vector(3*xstar_matrix-0.5)))

pd_T_dataframe_2 <- data.frame(matrix(NA, nrow = p, ncol = N))

for (r in 1:p) {
  pd_T_dataframe_2[r,] <-
    matrix(betahat_matrix[r+1,], nrow = N) + 
    (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (2*lhat_vector[r]^2)) *
       txstar_T_matrix) %*% 
    Ainv_matrix %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
  cat("pd_t:", r, "/", p, ",", sep = "")
}
points(as.numeric(as.vector(pd_T_dataframe_2[1,])) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "green")

A_matrix_3 <- matrix(NA, nrow = n, ncol = n)
for (r in 1:n) {
  for (c in 1:n) {
    if (is.na(A_matrix_3[r,c])) {
      A_matrix_3[r,c] <- A_matrix_3[c,r] <- 
        exp(-matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
              (lhatinv_diag_matrix^2) %*%
              matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
        )
    }
  }
}
Ainv_matrix_3 <- solve(A_matrix_3)

txstar_T_matrix_3 <- matrix(NA, nrow = N, ncol = n)
for (r in 1:N) {
  for (c in 1:n) {
    txstar_T_matrix_3[r,c] <- exp(-matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*% 
                                    (lhatinv_diag_matrix^2) %*%
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
    )
  }
  cat("txstarT: ", r, "/", N, ",", sep = "")
}

pd_T_dataframe_3 <- data.frame(matrix(NA, nrow = p, ncol = N))

for (r in 1:p) {
  pd_T_dataframe_3[r,] <-
    matrix(betahat_matrix[r+1,], nrow = N) + 
    2*(((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (lhat_vector[r]^2)) *
         txstar_T_matrix_3) %*% 
    Ainv_matrix_3 %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
  cat("pd_t:", r, "/", p, ",", sep = "")
}
points(as.numeric(as.vector(pd_T_dataframe_3[1,])) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "orange")

covMat1Mat2(object = )

pd_dataframe <- t(pd_T_dataframe)

pdnormed_dataframe <- pd_dataframe / sqrt(rowSums(pd_dataframe^2))

pdnormed_T_matrix <- t(data.matrix(pdnormed_dataframe))