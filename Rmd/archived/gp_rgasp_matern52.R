gp_rgasp_matern52 <- function(input_matrix, output_matrix, H_matrix){
  rgasp(design = input_matrix,
        response = output_matrix,
        trend = H_matrix,
        lower_bound=FALSE,
        zero.mean = "No"
  )
}