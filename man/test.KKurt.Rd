% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multNormTest.R
\name{test.KKurt}
\alias{test.KKurt}
\title{Test of normality based on Koziols measure of multivariate sample kurtosis}
\usage{
test.KKurt(data, MC.rep = 10000, alpha = 0.05)
}
\arguments{
\item{data}{a n x d matrix of d dimensional data vectors.}

\item{MC.rep}{number of repetitions for the Monte Carlo simulation of the critical value}

\item{alpha}{level of significance of the test}
}
\value{
a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
\describe{
        \item{\code{$Test}}{name of the test.}
        \item{\code{$Test.value}}{the value of the test statistic.}
        \item{\code{$cv}}{the approximated critical value.}
        \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.} }
}
\description{
Computes the multivariate normality test based on the invariant measure of multivariate sample kurtosis due to Koziol (1989).
}
\details{
Multivariate sample kurtosis due to Koziol (1989) is defined by
\deqn{\widetilde{b}_{n,d}^{(2)}=\frac{1}{n^2}\sum_{j,k=1}^n(Y_{n,j}^\top Y_{n,k})^4,}
where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
\eqn{n \ge d+1} is needed. If that is not the case the test returns an error. Note that for \eqn{d=1}, we have a measure proportional to the squared sample kurtosis.
}
\examples{
test.KKurt(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)

}
\references{
Koziol, J.A. (1989), A note on measures of multivariate kurtosis, Biom. J., 31:619-624.
}
\seealso{
\code{\link{KKurt}}
}
