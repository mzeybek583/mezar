find_peaks <- function(x, ...) {
  x <- as.vector(x)
  dens <- density(x, ...)
  
  second_deriv <- diff(sign(diff(dens$y)))
  dens$x[which(second_deriv == -2) + 1]
}