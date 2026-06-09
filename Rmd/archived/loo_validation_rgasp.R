loo_validation_rgasp <- function(cov_type){
  
  loo_df <- data.frame(matrix(NA, nrow = num_obs, ncol = 5))
  colnames(loo_df) <- c("actual", "loo_prediction", "lower95", "upper95", "validated")
  
  if(cov_type == "gauss"){
    for (obs in 1:num_obs) {
      gp_object <- gp_rgasp_gauss(input[-obs,,drop = F], output[-obs,,drop = F], H[-obs,,drop = F])
      loo_df$actual[obs] <- output[obs,,drop = F]
      loo_df$loo_prediction[obs] <- gppred_rgasp(gp_object, input[obs,,drop = F], H[obs,,drop = F])@mean
      loo_df$lower95[obs] <- gppred_rgasp(gp_object, input[obs,,drop = F], H[obs,,drop = F])@lower95
      loo_df$upper95[obs] <- gppred_rgasp(gp_object, input[obs,,drop = F], H[obs,,drop = F])@upper95
      loo_df$validated[obs] <- 
        loo_df$actual[obs] >= loo_df$lower95[obs] &
        loo_df$actual[obs] <= loo_df$upper95[obs]
    }
  }
  
  if(cov_type == "matern52"){
    for (obs in 1:num_obs) {
      gp_object <- gp_rgasp_matern52(input[-obs,,drop = F], output[-obs,,drop = F], H[-obs,,drop = F])
      loo_df$actual[obs] <- output[obs,,drop = F]
      loo_df$loo_prediction[obs] <- gppred_rgasp(gp_object, input[obs,,drop = F], H[obs,,drop = F])@mean
      loo_df$lower95[obs] <- gppred_rgasp(gp_object, input[obs,,drop = F], H[obs,,drop = F])@lower95
      loo_df$upper95[obs] <- gppred_rgasp(gp_object, input[obs,,drop = F], H[obs,,drop = F])@upper95
      loo_df$validated[obs] <- 
        loo_df$actual[obs] >= loo_df$lower95[obs] &
        loo_df$actual[obs] <= loo_df$upper95[obs]
    }
  }
  
  assign(paste0("gp_RG_", cov_type, "_loo"), 
         loo_df, 
         envir = .GlobalEnv)
}