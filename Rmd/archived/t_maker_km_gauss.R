t_maker_km_gauss <- function(inputs, test_inputs, range_vector){
  t_matrix <- matrix(NA, nrow = nrow(test_inputs), ncol = nrow(inputs))
  rangeinv_diag_matrix <- diag(x = 1/range_vector, 
                               nrow = length(range_vector)
  )
  for (r in 1:nrow(test_inputs)) {
    for (c in 1:nrow(inputs)) {
      t_matrix[r,c] <- exp(-0.5 *
                             matrix(test_inputs[r,] - inputs[c,], nrow = 1) %*% 
                                    (rangeinv_diag_matrix^2) %*%
                                    matrix(test_inputs[r,] - inputs[c,], ncol = 1)
      )
    }
    # cat("txstarT: ", r, "/", num_testing_input, ",", sep = "")
  }
  t_matrix
}