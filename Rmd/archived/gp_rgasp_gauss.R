gp_rgasp_gauss <- function(input_matrix, output_matrix, H_matrix){
  rgasp(design = input_matrix,
        response = output_matrix,
        trend = H_matrix,
        lower_bound=FALSE,
        kernel_type = 'pow_exp', alpha = rep(2, nrow(input_matrix)),
        zero.mean = "No"
        )
}