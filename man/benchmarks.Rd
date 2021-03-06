% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bencharks.R
\name{benchmarks}
\alias{benchmarks}
\title{Benchmark R and Rcpp functions.}
\description{
Use R package \code{microbenchmark} to compare the performance of C functions (\code{MetropolisR}) and Cpp functions (\code{MetropolisC}).
}
\examples{
\dontrun{
data(data)
attach(data)
tm <- microbenchmark::microbenchmark(
  rw <- rw.MetropolisR(.5, 25, 2000)
  rw <- rw.MetropolisC(.5, 25, 2000)
)
print(summary(tm)[,c(1,3,5,6)])
}
}
