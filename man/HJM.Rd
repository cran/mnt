% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{HJM}
\alias{HJM}
\title{statistic of the Henze-Jiménes-Gamero-Meintanis test}
\usage{
HJM(data, a)
}
\arguments{
\item{data}{a n x d numeric matrix of data values.}

\item{a}{positive numeric number (tuning parameter).}
}
\value{
The value of the test statistic.
}
\description{
Computes the test statistic of the Henze-Jiménes-Gamero-Meintanis test.
}
\details{
This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
\eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
}
\examples{
HJM(MASS::mvrnorm(20,c(0,1),diag(1,2)),a=2.5)

}
\references{
Henze, N., Jiménes-Gamero, M.D., Meintanis, S.G. (2019), Characterizations of multinormality and corresponding tests of fit, including for GARCH models, Econometric Th., 35:510–546, \href{https://doi.org/10.1017/S0266466618000154}{DOI}.
}
