gppred_rgasp <- function(rgasp_object, testing_inputs = testing_input, trend_matrix = H_testing) {
  RobustGaSP::predict(rgasp_object, testing_inputs, 
          testing_trend = trend_matrix,
          outasS3 = FALSE)
}