% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{DEHT}
\alias{DEHT}
\title{Statistic of the DEH test based on harmonic oscillator}
\usage{
DEHT(data, a = 1)
}
\arguments{
\item{data}{a n x d numeric matrix of data values.}

\item{a}{positive numeric number (tuning parameter).}
}
\value{
The value of the test statistic.
}
\description{
Computes the test statistic of the DEH test.
}
\details{
This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
\eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
}
\examples{
DEHT(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=1)

}
\references{
Dörr, P., Ebner, B., Henze, N. (2019) "Testing multivariate normality by zeros of the harmonic oscillator in characteristic function spaces" \href{https://arxiv.org/abs/1909.12624}{arXiv:1909.12624}
}
