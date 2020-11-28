
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

las <- segment_shapes(las, shp_plane(k = 10), "Coplanar")
plot(las, color = "Coplanar")

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

is.planar <- function(x, y, z, th1 = 12, th2 = 3) {
  xyz <- cbind(x,y,z)
  eigen_m <- eigen_values(xyz)
  is_planar <- eigen_m[2] > (th1*eigen_m[3]) && (th2*eigen_m[2]) > eigen_m[1]
  return(list(planar = is_planar))
}

M <- point_metrics(las, ~is.planar(X,Y,Z), k = 20)

las <- add_attribute(las, FALSE, "planar")
las$planar[M$pointID] <- M$planar
plot(las, color = "planar")
