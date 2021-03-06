% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{CS}
\alias{CS}
\title{Statistic of the test of Cox and Small}
\usage{
CS(data, Points = NULL)
}
\arguments{
\item{data}{a n x d matrix of d dimensional data vectors.}

\item{Points}{points for approximation of the maximum on the sphere. \code{Points=NULL} generates 5000 uniformly distributed Points on the d dimensional unit sphere.}
}
\value{
approximation of the value of the test statistic of the test of Cox and Small (1978).
}
\description{
This function returns the (approximated) value of the test statistic of the test of Cox and Small (1978).
}
\details{
The test statistic is \eqn{T_{n,CS}=\max_{b\in\{x\in\mathbf{R}^d:\|x\|=1\}}\eta_n^2(b)},
where \deqn{\eta_n^2(b)=\frac{\left\|n^{-1}\sum_{j=1}^nY_{n,j}(b^\top Y_{n,j})^2\right\|^2-\left(n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^3\right)^2}{n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^4-1-\left(n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^3\right)^2}}.
Here, \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
\eqn{n \ge d+1} is needed. If that is not the case the function returns an error. Note that the maximum functional has to be approximated by a discrete version, for details see Ebner (2012).
}
\examples{
CS(MASS::mvrnorm(50,c(0,1),diag(1,2)))

}
\references{
Cox, D.R. and Small, N.J.H. (1978), Testing multivariate normality, Biometrika, 65:263–272.

Ebner, B. (2012), Asymptotic theory for the test for multivariate normality by Cox and Small, Journal of Multivariate Analysis, 111:368–379.
}
