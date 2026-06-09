plot_preds_vs_byhand_DK <- function(inputs, outputs, 
                                    testing_inputs, prediction_object,
                                    byhand_predictions){
  # plot 1
  
  plot(NULL,
       xlim = c(start_point-amount,start_point+1+amount),
       ylim = c(-2,2),
       xlab="input",
       ylab="output"
  )
  legend("topleft",
         c("DiceKriging","by-hand","observations"),
         fill=c("black","red","blue")
  )
  #
  points(prediction_object$mean ~ testing_inputs, pch=20, cex=0.75)

  points(byhand_predictions ~ testing_inputs, pch=4, cex=0.25, col = "red")
  
  points(outputs ~ inputs, pch=19, col="blue", cex=1.2)
  
  # plot 2
  
  plot((prediction_object$mean - as.vector(byhand_predictions)) ~ testing_inputs,
       xlab = "input", ylab = "(DiceKriging prediction) - (by-hand prediction)")
  
  points(rep(0,num_obs) ~ inputs, pch = 19, col = "blue")
}