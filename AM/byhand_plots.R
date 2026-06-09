byhand_plots <- function(inputs, outputs, 
                                    testing_inputs, prediction_object,
                                    byhand_predictions){
  # plot 1
  
  plot(byhand_predictions ~ prediction_object$mean,
       col="black", pch=19, cex=0.25,
       xlab = "RobustGaSP predictions",
       ylab = "By-hand predictions")
  abline(0,1, col="blue")
  
  # plot 2
  
  plot((prediction_object$mean - as.vector(byhand_predictions)) ~ c(1:nrow(testing_inputs)),
       xlab = "test input index", ylab = "(RobustGaSP prediction) - (by-hand prediction)")
  abline(0,0, col="blue")
  
}