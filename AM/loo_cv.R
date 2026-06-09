loo_cv <- function(rgasp_object, outputs, gridbox_number){
  
  loo_df <- data.frame(matrix(NA, nrow = num_obs, ncol = 5))
  colnames(loo_df) <- c("actual", "loo_prediction", "lower95", "upper95", "validated")
  
  loo_df$actual <- outputs
  loo_df$loo_prediction <- leave_one_out_rgasp(rgasp_object)$mean
  loo_df$lower95 <- leave_one_out_rgasp(rgasp_object)$mean - qnorm(0.975) * leave_one_out_rgasp(rgasp_object)$sd
  loo_df$upper95 <- leave_one_out_rgasp(rgasp_object)$mean + qnorm(0.975) * leave_one_out_rgasp(rgasp_object)$sd
  loo_df$validated <- loo_df$actual >= loo_df$lower95 & loo_df$actual <= loo_df$upper95
  
  assign(paste0("loo_cv", gridbox_number), 
         loo_df, 
         envir = .GlobalEnv)
}