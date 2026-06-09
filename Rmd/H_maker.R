H_maker <- function(inputs) {
  cbind(c(rep(1, nrow(inputs))), inputs)
}