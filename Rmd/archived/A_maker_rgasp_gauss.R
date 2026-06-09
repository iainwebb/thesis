A_maker_rgasp_gauss <- function(inputs, range_vector){
  A_matrix <- matrix(NA, nrow = nrow(inputs), ncol = nrow(inputs))
  rangeinv_diag_matrix <- diag(x = 1/range_vector, 
                               nrow = length(range_vector)
                               )
  for (r in 1:nrow(inputs)) {
    for (c in 1:nrow(inputs)) {
      if (is.na(A_matrix[r,c])) {
        A_matrix[r,c] <- A_matrix[c,r] <- 
          exp(-1*
                matrix(inputs[r,] - inputs[c,], nrow = 1) %*%
                (rangeinv_diag_matrix^2) %*%
                matrix(inputs[r,] - inputs[c,], ncol = 1)
              )
      }
    }
  }
  A_matrix
}