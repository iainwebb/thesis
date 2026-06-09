m_map_loo <- function(inputs, outputs, cov_type){
  
  loo_df <- data.frame(matrix(NA, nrow = num_obs, ncol = 5))
  colnames(loo_df) <- c("actual", "loo_prediction", "lower95", "upper95", "validated")
  
  if(cov_type == "gauss"){
    for (obs in 1:num_obs) {
      gp_object <- rgasp(design = inputs[-obs,,drop = F], response = outputs[-obs,,drop = F], trend = H[-obs,,drop = F], 
                         lower_bound=FALSE, zero.mean = "No",
                         kernel_type = 'pow_exp', alpha = rep(2, nrow(inputs)),        # squared exponential
      )
      loo_df$actual[obs] <- outputs[obs,,drop = F]
      loo_df$loo_prediction[obs] <- gppred_rgasp(gp_object, inputs[obs,,drop = F], H[obs,,drop = F])@mean
      loo_df$lower95[obs] <- gppred_rgasp(gp_object, inputs[obs,,drop = F], H[obs,,drop = F])@lower95
      loo_df$upper95[obs] <- gppred_rgasp(gp_object, inputs[obs,,drop = F], H[obs,,drop = F])@upper95
      loo_df$validated[obs] <- 
        loo_df$actual[obs] >= loo_df$lower95[obs] &
        loo_df$actual[obs] <= loo_df$upper95[obs]
    }
  }
  
  if(cov_type == "matern52"){
    for (obs in 1:num_obs) {
      gp_object <- rgasp(design = inputs[-obs,,drop = F], response = outputs[-obs,,drop = F], trend = H[-obs,,drop = F], 
                         lower_bound=FALSE, zero.mean = "No"
      )
      loo_df$actual[obs] <- outputs[obs,,drop = F]
      loo_df$loo_prediction[obs] <- gppred_rgasp(gp_object, inputs[obs,,drop = F], H[obs,,drop = F])@mean
      loo_df$lower95[obs] <- gppred_rgasp(gp_object, inputs[obs,,drop = F], H[obs,,drop = F])@lower95
      loo_df$upper95[obs] <- gppred_rgasp(gp_object, inputs[obs,,drop = F], H[obs,,drop = F])@upper95
      loo_df$validated[obs] <- 
        loo_df$actual[obs] >= loo_df$lower95[obs] &
        loo_df$actual[obs] <= loo_df$upper95[obs]
    }
  }
  
  assign(paste0("m.map.", cov_type, ".loo"), 
         loo_df, 
         envir = .GlobalEnv)
}