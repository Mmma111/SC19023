#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp
//' @param sigma The standard deviation of proposal distribution N(Xt,sigma^2)
//' @param X0  initial value X0
//' @param N the length of the chain
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rwC <- rw.MetropolisC(.5, 25, 2000)
//'  plot(a, rwC$x, type = "s", main = "", xlab = "sigma = 0.5", ylab = "x")
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolisC(double sigma, double x0, double N) {
  NumericVector x(N+1);
  x[0] = x0;
  NumericVector u = runif(N);
  double k = 0;
  double y = 0;
  for (int i = 2; i < N+1; i++) {
    y = rnorm(1, x[i-2], sigma)[0];
    if (u[i-2] <= exp(-((abs(y))-(abs(x[i-2])))) ) {
      x[i-1] = y;
    } 
    else {
      x[i-1] = x[i-2];
      k++;
    }
  }
  x[N] = k;
  return(x);
}

