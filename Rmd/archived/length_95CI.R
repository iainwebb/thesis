length_95CI <- function(posterior_lower95, posterior_upper95){
  sum(posterior_upper95 - posterior_lower95)/num_testing_input
}