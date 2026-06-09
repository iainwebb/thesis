loo_validation_plot <- function(loo_df){
  num_obs <- nrow(loo_df)
  cols <- rep("red", num_obs)
  cols[loo_df$validated] <- "black"
  
  plot(loo_df$loo_prediction ~ loo_df$actual,
       xlim = c(min(loo_df$lower95, loo_df$actual), 
                max(loo_df$upper95, loo_df$actual)),
       ylim = c(min(loo_df$lower95, loo_df$actual), 
                max(loo_df$upper95, loo_df$actual)),
       col = cols, pch = 20,
       xlab = "Model", ylab = "Emulator")
  abline(0,1, col="blue")
  segments(loo_df$actual,
           loo_df$lower95,
           loo_df$actual,
           loo_df$upper95,
           col = cols)
}