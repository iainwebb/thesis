byhand_preds <- function(testing_inputs, 
                            trend_vector, 
                            t_matrix, 
                            A_inv_matrix,
                            outputs,
                            H_matrix){
  matrix(
    c(rep(1, nrow(testing_inputs)), testing_inputs), 
    nrow = nrow(testing_inputs), byrow = F
    ) %*% trend_vector +
    t_matrix %*% 
    A_inv_matrix %*% 
    (as.matrix(outputs, ncol = nrow(outputs)) - H_matrix %*% trend_vector)
}