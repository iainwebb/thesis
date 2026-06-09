gp_km_matern52 <- function(input_matrix, output_matrix){
  km(~.,
     input_matrix,
     output_matrix,
     covtype="matern5_2",
     optim.method="BFGS",
     control=list(maxit=500)
     )
}