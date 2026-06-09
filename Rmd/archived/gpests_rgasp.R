gpests_rgasp <- function(rgasp_object, cov_type){
  ifelse(cov_type == "gauss",
         assign(paste0("gp_rg_", cov_type, "_range_pars"), 
                (1 / rgasp_object@beta_hat), 
                envir = .GlobalEnv),
         assign(paste0("gp_rg_", cov_type, "_range_pars"), 
                1 / rgasp_object@beta_hat, 
                envir = .GlobalEnv)
         )
  assign(paste0("gp_rg_", cov_type, "_sd_par"), 
         rgasp_object@sigma2_hat, 
         envir = .GlobalEnv)
  assign(paste0("gp_rg_", cov_type, "_trend_pars"), 
         rgasp_object@theta_hat,
         envir = .GlobalEnv)
}