#------------------------
# a 3 dimensional example
#------------------------
# dimensional of the inputs
dim_inputs <- 3
# number of the inputs
num_obs <- 30
# uniform samples of design
# input <- matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs)
# Following codes use maximin Latin Hypercube Design, which is typically better than uniform
library(lhs)
input <- maximinLHS(n=num_obs, k=dim_inputs) ##maximin lhd sample
####
library(RobustGaSP)
# outputs from the 3 dim dettepepel.3.data function
output = matrix(0,num_obs,1)
for(i in 1:num_obs){
  output[i]<-dettepepel.3.data(input[i,])
}
# use constant mean basis, with no constraint on optimization
m1<- rgasp(design = input, response = output, lower_bound=FALSE)
# the following use constraints on optimization
# m1<- rgasp(design = input, response = output, lower_bound=TRUE)
# the following use a single start on optimization
# m1<- rgasp(design = input, response = output, lower_bound=FALSE)
# number of points to be predicted
num_testing_input <- 5000
# generate points to be predicted
testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
# Perform prediction
m1.predict<-predict(m1, testing_input, outasS3 = FALSE)
# Predictive mean
#m1.predict@mean
# The following tests how good the prediction is
testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-dettepepel.3.data(testing_input[i,])
}
# compute the MSE, average coverage and average length
# out of sample MSE
MSE_emulator <- sum((m1.predict@mean-testing_output)^2)/(num_testing_input)
# proportion covered by 95% posterior predictive credible interval
prop_emulator <- length(which((m1.predict@lower95<=testing_output)
                              &(m1.predict@upper95>=testing_output)))/num_testing_input
# average length of posterior predictive credible interval
length_emulator <- sum(m1.predict@upper95-m1.predict@lower95)/num_testing_input
# output of prediction
MSE_emulator
prop_emulator
length_emulator
# normalized RMSE
sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))

## 1D example
n <- 5
inputs_x_norm_T_matrix <- randomLHS(n,1)
inputs_x_norm_matrix <- t(inputs_x_norm_T_matrix)

## There are p variables in all, and n runs
p <- ncol(inputs_x_norm_T_matrix)
n <- nrow(inputs_x_norm_T_matrix)

response_test_dataframe <- 
  data.frame(y = matrix(2*sin(3*inputs_x_norm_T_matrix - 0.5), nrow = n))
plot(response_test_dataframe$y ~ as.vector(inputs_x_norm_matrix))
# GP set-up ####################################################################

## We specify some variables, namely N (the size of the sample from the input
## parameter space) and q (the number of betas that will be estimated) 

N <- 10
q <- 1 + p

## We load the required packages for sampling
library(lhs)

## We create the sample
xstar_matrix <- t(randomLHS(N, p))
# TO DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# xstar_T_dataframe <- data.frame(t(xstar_matrix))
# colnames(xstar_T_dataframe) <- colnames(inputs_x_norm_T_matrix)
# TO DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Our estimated posterior derivative requires H, which doesn't depend on the 
# response (and thus on the gridbox in the case of AOD)
H_matrix <- cbind(c(rep(1, n)), t(unname(inputs_x_norm_matrix)))

# use constant mean basis, with no constraint on optimization
library(RobustGaSP)
m1<- rgasp(design = data.frame(inputs_x_norm_T_matrix), 
           response = response_test_dataframe, 
           lower_bound=FALSE
           )
# the following use constraints on optimization
# m1<- rgasp(design = input, response = output, lower_bound=TRUE)
# the following use a single start on optimization
# m1<- rgasp(design = input, response = output, lower_bound=FALSE)
# number of points to be predicted
num_testing_input <- 5000
# generate points to be predicted
testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
# Perform prediction
m1.predict<-predict(m1, testing_input, outasS3 = FALSE)
# Predictive mean
#m1.predict@mean
# The following tests how good the prediction is
testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-dettepepel.3.data(testing_input[i,])
}
# compute the MSE, average coverage and average length
# out of sample MSE
MSE_emulator <- sum((m1.predict@mean-testing_output)^2)/(num_testing_input)
# proportion covered by 95% posterior predictive credible interval
prop_emulator <- length(which((m1.predict@lower95<=testing_output)
                              &(m1.predict@upper95>=testing_output)))/num_testing_input
# average length of posterior predictive credible interval
length_emulator <- sum(m1.predict@upper95-m1.predict@lower95)/num_testing_input
# output of prediction
MSE_emulator
prop_emulator
length_emulator
# normalized RMSE
sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))
