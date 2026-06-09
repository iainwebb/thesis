gp_km_gauss <- function(input_matrix, output_matrix){
  km(~.,
     input_matrix,
     output_matrix,
     covtype="gauss",
     optim.method="BFGS",
     control=list(maxit=500)
     )
}