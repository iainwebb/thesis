loo_validation_km <- function(inputs, outputs, cov_type){
  
  num_obs <- nrow(inputs)
  
  loo_df <- data.frame(matrix(NA, nrow = num_obs, ncol = 5))
  colnames(loo_df) <- c("actual", "loo_prediction", "lower95", "upper95", "validated")
  
  if(cov_type == "gauss"){
    for (obs in 1:num_obs) {
      gp_object <- gp_km_gauss(inputs[-obs,,drop = F], outputs[-obs,,drop = F])
      loo_df$actual[obs] <- outputs[obs,,drop = F]
      loo_df$loo_prediction[obs] <- gppred_km(gp_object, inputs[obs,,drop = F])$mean
      loo_df$lower95[obs] <- gppred_km(gp_object, inputs[obs,,drop = F])$lower95
      loo_df$upper95[obs] <- gppred_km(gp_object, inputs[obs,,drop = F])$upper95
      loo_df$validated[obs] <- 
        loo_df$actual[obs] >= loo_df$lower95[obs] &
        loo_df$actual[obs] <= loo_df$upper95[obs]
    }
  }
  
  if(cov_type == "matern52"){
    for (obs in 1:num_obs) {
      gp_object <- gp_km_matern52(inputs[-obs,,drop = F], outputs[-obs,,drop = F])
      loo_df$actual[obs] <- outputs[obs,,drop = F]
      loo_df$loo_prediction[obs] <- gppred_km(gp_object, inputs[obs,,drop = F])$mean
      loo_df$lower95[obs] <- gppred_km(gp_object, inputs[obs,,drop = F])$lower95
      loo_df$upper95[obs] <- gppred_km(gp_object, inputs[obs,,drop = F])$upper95
      loo_df$validated[obs] <- 
        loo_df$actual[obs] >= loo_df$lower95[obs] &
        loo_df$actual[obs] <= loo_df$upper95[obs]
    }
  }

  assign(paste0("gp_DK_", cov_type, "_loo"), 
         loo_df, 
         envir = .GlobalEnv)
}