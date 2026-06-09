gppred_km <- function(km_object, testing_inputs = testing_input) {
  predict(km_object, testing_inputs, "SK", checkNames = F)
}