estimated_parameters <- function(emulator, gridbox_number){
  
  assign(paste0("estimated_trends", gridbox_number), 
         emulator@theta_hat,
         envir = .GlobalEnv)
  assign(paste0("estimated_ranges", gridbox_number), 
         (1 / emulator@beta_hat), 
         envir = .GlobalEnv)
  assign(paste0("estimated_var", gridbox_number), 
         emulator@sigma2_hat, 
         envir = .GlobalEnv)
}