#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{MetropolisR}) and Cpp functions (\code{MetropolisC}).
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma
#' @useDynLib SC19023
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' tm <- microbenchmark::microbenchmark(
#'   rw <- rw.MetropolisR(.5, 25, 2000)
#'   rw <- rw.MetropolisC(.5, 25, 2000)
#' )
#' print(summary(tm)[,c(1,3,5,6)])
#' }
NULL