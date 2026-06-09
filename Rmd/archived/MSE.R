MSE <- function(posterior_mean, actual_output){
  sum((posterior_mean-actual_output)^2)/(num_testing_input)
}