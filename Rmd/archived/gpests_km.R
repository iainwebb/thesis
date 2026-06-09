gpests_km <- function(km_object, cov_type){
  assign(paste0("gp_dk_", cov_type, "_range_pars"), 
         km_object@covariance@range.val, 
         envir = .GlobalEnv)
  assign(paste0("gp_dk_", cov_type, "_sd_par"), 
         km_object@covariance@sd2, 
         envir = .GlobalEnv)
  assign(paste0("gp_dk_", cov_type, "_trend_pars"), 
         km_object@trend.coef,
         envir = .GlobalEnv)
}