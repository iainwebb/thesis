## We define a function called gp_d_dxi_normalised_function which returns the
## normalised GP estimates of the partial derivatives. It takes as its single
## argument a 221-long vector or 221 x 1 dataframe of responses

gp_d_dxi_normalised_function <- function(response, loglik_index) {
  
  gp <- km(~., 
           design = matrix(inputs_x_norm_T_matrix), 
           response = response_test_dataframe, 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
  )
  
  ## Using the posterior estimates of the range parameters, we make the matrix 
  ## A, its inverse, and t(xstar)^T
  lhat_vector <- gp@covariance@range.val
  lhat_matrix <- matrix(lhat_vector)
  lhatinv_diag_matrix <- diag(x = 1/lhat_vector, nrow = p)
  sqrt(gp@covariance@sd2)
  sqrt(2)
  A_matrix <- matrix(NA, nrow = n, ncol = n)
  for (r in 1:n) {
    for (c in 1:n) {
      if (is.na(A_matrix[r,c])) {
        A_matrix[r,c] <- A_matrix[c,r] <- 
          exp(-0.5*
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                (lhatinv_diag_matrix^2) %*%
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
          )
      }
    }
  }
  Ainv_matrix <- solve(A_matrix)
  
  txstar_T_matrix <- matrix(NA, nrow = N, ncol = n)
  for (r in 1:N) {
    for (c in 1:n) {
      txstar_T_matrix[r,c] <- exp(-0.5 *
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*% 
                                    (lhatinv_diag_matrix^2) %*%
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
      )
    }
    cat("txstarT: ", r, "/", N, ",", sep = "")
  }
  
  betahat_matrix <- matrix(gp@trend.coef, nrow = q)
  
  plot(response_test_dataframe$y ~ as.vector(inputs_x_norm_matrix), xlim = c(0,1), ylim = c(-2,2))
  
  predictions <- predict(gp, as.vector(t(xstar_matrix)), "SK")
  points(predictions$mean ~ as.vector(xstar_matrix), pch=19, cex=0.5)
  
  my_predictions <- matrix(c(rep(1, N), xstar_matrix), nrow = N, byrow = F) %*% betahat_matrix +
                      txstar_T_matrix %*% Ainv_matrix %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
  points(my_predictions ~ as.vector(xstar_matrix[1,]), pch=19, cex=0.5, col = "red")
  
  sum(abs(predictions$mean - as.vector(my_predictions))) / (10^(-8)*N)
  plot((predictions$mean - as.vector(my_predictions)) ~ as.vector(xstar_matrix))
  points(rep(0,5) ~ as.vector(inputs_x_norm_matrix), pch = 19, col = "red")
  
  ## Secondly, we use the trend estimates to estimate the 37 posterior partial
  ## derivatives at each of the N points, then normalise them
  
  pd_T_dataframe <- data.frame(matrix(NA, nrow = p, ncol = N))
  
  for (r in 1:p) {
    pd_T_dataframe[r,] <-
      matrix(betahat_matrix[r+1,], nrow = N) + 
      (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (lhat_vector[r]^2)) *
         txstar_T_matrix) %*% 
      Ainv_matrix %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
    cat("pd_t:", r, "/", p, ",", sep = "")
  }


  plot(response_test_dataframe$y ~ as.vector(inputs_x_norm_matrix), xlim = c(0,1), ylim = c(-20,20))
  
  points((6*cos(as.vector(3*xstar_matrix-0.5))) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "blue")
  points(as.numeric(as.vector(pd_T_dataframe[1,])) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "red")

  as.numeric(as.vector(pd_T_dataframe[1,])) / (6*cos(as.vector(3*xstar_matrix-0.5)))

  pd_T_dataframe_2 <- data.frame(matrix(NA, nrow = p, ncol = N))
  
  for (r in 1:p) {
    pd_T_dataframe_2[r,] <-
      matrix(betahat_matrix[r+1,], nrow = N) + 
      (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (2*lhat_vector[r]^2)) *
         txstar_T_matrix) %*% 
      Ainv_matrix %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
    cat("pd_t:", r, "/", p, ",", sep = "")
  }
  points(as.numeric(as.vector(pd_T_dataframe_2[1,])) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "green")
  
  A_matrix_3 <- matrix(NA, nrow = n, ncol = n)
  for (r in 1:n) {
    for (c in 1:n) {
      if (is.na(A_matrix_3[r,c])) {
        A_matrix_3[r,c] <- A_matrix_3[c,r] <- 
          exp(-matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                (lhatinv_diag_matrix^2) %*%
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
          )
      }
    }
  }
  Ainv_matrix_3 <- solve(A_matrix_3)
  
  txstar_T_matrix_3 <- matrix(NA, nrow = N, ncol = n)
  for (r in 1:N) {
    for (c in 1:n) {
      txstar_T_matrix_3[r,c] <- exp(-matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*% 
                                    (lhatinv_diag_matrix^2) %*%
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
      )
    }
    cat("txstarT: ", r, "/", N, ",", sep = "")
  }
  
  pd_T_dataframe_3 <- data.frame(matrix(NA, nrow = p, ncol = N))
  
  for (r in 1:p) {
    pd_T_dataframe_3[r,] <-
      matrix(betahat_matrix[r+1,], nrow = N) + 
      2*(((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (lhat_vector[r]^2)) *
         txstar_T_matrix_3) %*% 
      Ainv_matrix_3 %*% (as.matrix(response_test_dataframe, ncol = n) - H_matrix %*% betahat_matrix)
    cat("pd_t:", r, "/", p, ",", sep = "")
  }
  points(as.numeric(as.vector(pd_T_dataframe_3[1,])) ~ as.vector(xstar_matrix), pch=19, cex=0.5, col = "orange")
  
  covMat1Mat2(object = )
  
  pd_dataframe <- t(pd_T_dataframe)
  
  pdnormed_dataframe <- pd_dataframe / sqrt(rowSums(pd_dataframe^2))
  
  pdnormed_T_matrix <- t(data.matrix(pdnormed_dataframe))
  
  ## just when repeating runs
  # assign(paste0("loglik_long142lat88_", loglik_index), logLik.km(gp), 
  #        envir = .GlobalEnv)
  
  return(pdnormed_T_matrix)
}