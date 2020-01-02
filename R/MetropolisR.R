#' @title A random walk Metropolis sampler for generating the standard Laplace distribution using R
#' @description A random walk Metropolis sampler for generating the standard Laplace distribution using R
#' @importFrom stats runif rnorm
#' @param sigma The standard deviation of proposal distribution N(Xt,sigma^2)
#' @param X0  initial value X0
#' @param N the length of the chain
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' rwR <- rw.MetropolisR(.5, 25, 2000)
#' plot(a, rwR$x, type = "s", main = "", xlab = "sigma = 0.5", ylab = "x")
#' }
#' @export
rwMetropolisR <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(abs(x[i-1]) - abs(y)))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x))
}