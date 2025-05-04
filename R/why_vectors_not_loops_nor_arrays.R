# basket <- data.frame("Food" = c("Apples", "Bananas", "Carrots"),
                     # "PricePerUnit" = c(.99, .19, .49),
                     # "Quantity" = c(12, 6, 2))
basket <- rbind(basket, basket)
nrow(basket)


# Using a for loop ---------
# create the total
total <- 0
# loop over the data.frame and add the running total
start_time <- Sys.time()
for (row in 1:nrow(basket)) {
  total <- total + (basket$PricePerUnit[row] * basket$Quantity[row])
}
total
end_time <- Sys.time()
end_time - start_time

# Using apply --------------
# define the multiplication function
# multiply each PricePerUnit and Quantity and store the resulting vector
bsktTotal <- function(bskt){
  as.numeric(unlist((bskt)["PricePerUnit"])) * as.numeric(unlist((bskt)["Quantity"]))
}
perItemTotal <- apply(basket, 1, bsktTotal)
# sum all values in the perItemTotal
sum(perItemTotal)

# multiply each PricePerUnit and Quantity and store the resulting vector
start_time <- Sys.time()
perItemTotal <- apply(basket, 1, function(bskt) {
  as.numeric(unlist((bskt)["PricePerUnit"])) * as.numeric(unlist((bskt)["Quantity"]))
}
)
# sum all values in the perItemTotal
sum(perItemTotal)
end_time <- Sys.time()
end_time - start_time

# Vector operations --------
# take the sum of multiplying PerPriceUnit and Quantity to get total cost
start_time <- Sys.time()
sum(basket$PricePerUnit * basket$Quantity)
end_time <- Sys.time()
end_time - start_time
