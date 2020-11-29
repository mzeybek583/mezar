
#

## Linearity calculation

# Library ---------------------------------------------------------------

library(lidR)

# Read data ---------------------------------------------------------------

data <- readLAS(files = "data/tek6-s3.las")

plot(data)

las = voxelize_points(data, 0.03)

las <- segment_shapes(las, algorithm = shp_line(th1=2, k = 50))

length(las@data$X[las@data$Shape==TRUE])

plot(las, color = "Shape", legend=TRUE)

Rcpp::sourceCpp(code = "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec eigen_values(arma::mat A) {
arma::mat coeff, score;
arma::vec latent;
arma::princomp(coeff, score, latent, A);
return(latent);
}")

is.linear <- function(x, y, z) {
  xyz <- cbind(x,y,z)
  eigen_m <- eigen_values(xyz)
  is_linear <- eigen_m[1] - eigen_m[2] / eigen_m[1]
  return(list(linear = is_linear))
}

M <- point_metrics(las, ~is.linear(X,Y,Z), k = 50)

las <- add_attribute(las, FALSE, "linear")
las$linear<- M$linear
plot(las, color = "linear", legend=T)
