prop_wn_95CI <- function(posterior_lower95, posterior_upper95){
  length(
    which((posterior_lower95<=testing_output)
               &(posterior_upper95>=testing_output))
    )/num_testing_input
}