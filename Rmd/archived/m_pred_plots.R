m_pred_plots <- function(inputs, outputs, 
                                    testing_inputs, prediction_object,
                                    byhand_predictions){
  # plot 1
  
  plot(NULL,
       xlim = c(start_point-amount,start_point+1+amount),
       ylim = c(min(prediction_object$mean),max(prediction_object$mean)),
       xlab="input",
       ylab="output"
  )
  legend("topleft",
         c("RobustGaSP","by-hand","observations"),
         fill=c("black","red","blue")
  )
  #
  points(prediction_object$mean ~ testing_input, pch=20, cex=0.75)

  points(byhand_predictions ~ testing_input, pch=4, cex=0.25, col = "red")
  
  points(outputs ~ inputs, pch=19, col="blue", cex=1.2)
  
  # plot 2
  
  plot((prediction_object$mean - as.vector(byhand_predictions)) ~ testing_inputs,
       xlab = "input", ylab = "(RobustGaSP prediction) - (by-hand prediction)")
  
  points(rep(0,num_obs) ~ inputs, pch = 19, col = "blue")
}