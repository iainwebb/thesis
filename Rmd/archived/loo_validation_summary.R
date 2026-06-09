loo_validation_summary <- function(loo_df, estimation_type, cov_type){
  kable(round(loo_df,2), 
        align=c("c","c","c","c"),
        caption = paste0(estimation_type, " with ", cov_type, " kernel (", 100*sum(loo_df[,5]) / num_obs, "% validation)")
  )
}