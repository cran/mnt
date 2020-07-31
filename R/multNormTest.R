# Code for concrete tests

#' Test of normality based on Mardias measure of multivariate sample skewness
#' @description
#' Computes the multivariate normality test based on the classical invariant measure of multivariate sample skewness due to Mardia (1970).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.} }
#'
#'
#' @details
#' Multivariate sample skewness due to Mardia (1970) is defined by
#' \deqn{b_{n,d}^{(1)}=\frac{1}{n^2}\sum_{j,k=1}^n(Y_{n,j}^\top Y_{n,k})^3,}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error. Note that for \eqn{d=1}, we have a measure proportional to the squared sample skewness.
#'
#' @references
#' Mardia, K.V. (1970), Measures of multivariate skewness and kurtosis with applications, Biometrika, 57:519-530.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467-506.
#'
#' @examples
#' test.MSkew(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{MSkew}}
#'
#' @export
test.MSkew<-function(data,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MSkew,tuning=NULL,repetitions = MC.rep)
  Testst=MSkew(data)
  result <- list("Test" = "MSkew", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}


#' Test of normality based on Mardias measure of multivariate sample kurtosis
#' @description
#' Computes the multivariate normality test based on the classical invariant measure of multivariate sample kurtosis due to Mardia (1970).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.} }
#'
#'
#' @details
#' Multivariate sample kurtosis due to Mardia (1970) is defined by
#' \deqn{b_{n,d}^{(2)}=\frac{1}{n}\sum_{j=1}^n\|Y_{n,j}\|^4,}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}.To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @references
#' Mardia, K.V. (1970), Measures of multivariate skewness and kurtosis with applications, Biometrika, 57:519-530.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467-506.
#'
#' @examples
#' test.MKurt(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{MKurt}}
#'
#' @export
test.MKurt<-function(data,MC.rep=10000,alpha=0.05){
    n=dim(data)[1]
    d=dim(data)[2]
    if (is.null(n)){if (is.vector(data)){
      n=length(data)
      d=1
    }else{warning("Wrong dimensions of data!")}}
    cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MKurt,tuning=NULL,repetitions = MC.rep)
    Testst=MKurt(data)
    result <- list("Test" = "MKurt", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
    attr(result, "class") <- "mnt"
    return(result)
  }



#' Test of normality based on Koziols measure of multivariate sample kurtosis
#' @description
#' Computes the multivariate normality test based on the invariant measure of multivariate sample kurtosis due to Koziol (1989).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.} }
#'
#'
#' @details
#' Multivariate sample kurtosis due to Koziol (1989) is defined by
#' \deqn{\widetilde{b}_{n,d}^{(2)}=\frac{1}{n^2}\sum_{j,k=1}^n(Y_{n,j}^\top Y_{n,k})^4,}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error. Note that for \eqn{d=1}, we have a measure proportional to the squared sample kurtosis.
#'
#' @references
#' Koziol, J.A. (1989), A note on measures of multivariate kurtosis, Biom. J., 31:619-624.
#'
#' @examples
#' test.KKurt(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{KKurt}}
#'
#' @export
test.KKurt<-function(data,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=KKurt,tuning=NULL,repetitions = MC.rep)
  Testst=KKurt(data)
  result <- list("Test" = "KKurt", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}





#' Test of normality based on multivariate skewness in the sense of Malkovich and Afifi
#' @description
#' Computes the test of multivariate normality based on skewness in the sense of Malkovich and Afifi (1973).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#' @param num.points number of points distributed uniformly over the sphere for approximation of the maximum on the sphere.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{number of points used in approximation.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#' }
#'
#'
#' @details
#' Multivariate sample skewness due to Malkovich and Afifi (1973) is defined by
#' \deqn{b_{n,d,M}^{(1)}=\max_{u\in \{x\in\mathbf{R}^d:\|x\|=1\}}\frac{\left(\frac{1}{n}\sum_{j=1}^n(u^\top X_j-u^\top \overline{X}_n )^3\right)^2}{(u^\top S_n u)^3},}
#' where \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @references
#' Malkovich, J.F., and Afifi, A.A. (1973), On tests for multivariate normality, J. Amer. Statist. Ass., 68:176-179.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467-506.
#'
#' @examples
#' \donttest{test.MASkew(MASS::mvrnorm(10,c(0,1),diag(1,2)),MC.rep=100)}
#'
#' @seealso
#' \code{\link{MASkew}}
#'
#' @export
test.MASkew<-function(data,MC.rep=10000,alpha=0.05,num.points = 1000){
    n=dim(data)[1]
    d=dim(data)[2]
    if (is.null(n)){if (is.vector(data)){
      n=length(data)
      d=1
    }else{warning("Wrong dimensions of data!")}}
    cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MASkew,tuning=GVP(num.points,d),repetitions = MC.rep)
    Testst=MASkew(data,GVP(num.points,d))
    result <- list("Test" = "MASkew", "param" = num.points, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
    attr(result, "class") <- "mnt"
    return(result)
  }




#' Test of normality based on multivariate kurtosis in the sense of Malkovich and Afifi
#' @description
#' Computes the multivariate normality test based on the invariant measure of multivariate sample kurtosis due to Malkovich and Afifi (1973).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#' @param num.points number of points distributed uniformly over the sphere for approximation of the maximum on the sphere.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{number of points used in approximation.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#'
#' @details
#' Multivariate sample skewness due to Malkovich and Afifi (1973) is defined by
#' \deqn{b_{n,d,M}^{(1)}=\max_{u\in \{x\in\mathbf{R}^d:\|x\|=1\}}\frac{\left(\frac{1}{n}\sum_{j=1}^n(u^\top X_j-u^\top \overline{X}_n )^3\right)^2}{(u^\top S_n u)^3},}
#' where \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @references
#' Malkovich, J.F., and Afifi, A.A. (1973), On tests for multivariate normality, J. Amer. Statist. Ass., 68:176-179.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467-506.
#'
#' @examples
#' \donttest{test.MAKurt(MASS::mvrnorm(10,c(0,1),diag(1,2)),MC.rep=100)}
#'
#' @seealso
#' \code{\link{MAKurt}}
#'
#' @export
test.MAKurt<-function(data,MC.rep=10000,alpha=0.05,num.points = 1000){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MAKurt,tuning=GVP(num.points,d),repetitions = MC.rep)
  Testst=MAKurt(data,GVP(num.points,d))
  result <- list("Test" = "MAKurt", "param" = num.points, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}


#' Test of multivariate normality based on the measure of multivariate skewness of Mori, Rohatgi and Szekely
#' @description
#' Computes the multivariate normality test based on the invariant measure of multivariate sample skewness due to Mori, Rohatgi and Szekely (1993).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.} }
#'
#'
#' @details
#' Multivariate sample skewness due to Mori, Rohatgi and Szekely (1993) is defined by
#' \deqn{\widetilde{b}_{n,d}^{(1)}=\frac{1}{n}\sum_{j=1}^n\|Y_{n,j}\|^2\|Y_{n,k}\|^2Y_{n,j}^\top Y_{n,k},}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error. Note that for \eqn{d=1}, it is equivalent to skewness in the sense of Mardia.
#'
#' @references
#' Mori, T. F., Rohatgi, V. K., Szekely, G. J. (1993), On multivariate skewness and kurtosis, Theory of Probability and its Applications, 38:547-551.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467-506.
#'
#' @examples
#' test.MRSSkew(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{MRSSkew}}
#'
#' @export
test.MRSSkew<-function(data,MC.rep=10000,alpha=0.05)
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MRSSkew,tuning=NULL,repetitions = MC.rep)
  Testst=MRSSkew(data)
  result <- list("Test" = "MRSSkew", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Baringhaus-Henze-Epps-Pulley (BHEP) test
#' @description
#' Performs the BHEP test of multivariate normality as suggested in Henze and Wagner (1997) using a tuning parameter \code{a}.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#'
#' @details
#' The test statistic is \deqn{BHEP_{n,\beta}=\frac{1}{n} \sum_{j,k=1}^n \exp\left(-\frac{\beta^2\|Y_{n,j}-Y_{n,k}\|^2}{2}\right)- \frac{2}{(1+\beta^2)^{d/2}} \sum_{j=1}^n \exp\left(- \frac{\beta^2\|Y_{n,j}\|^2}{2(1+\beta^2)} \right) + \frac{n}{(1+2\beta^2)^{d/2}}.}
#' Here, \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @references
#' Henze, N., Wagner, T. (1997), A new approach to the class of BHEP tests for multivariate normality, J. Multiv. Anal., 62:1-23, \href{https://doi.org/10.1006/jmva.1997.1684}{DOI}
#'
#' @examples
#' test.BHEP(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{BHEP}}
#'
#' @export
test.BHEP<-function(data,a=1,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=BHEP,tuning=a,repetitions = MC.rep)
  Testst=BHEP(data,a)
  result <- list("Test" = "BHEP", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' The Henze-Zirkler test
#' @description
#' Performs the test of multivariate normality of Henze and Zirkler (1990).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#'
#' @details
#' A \code{\link{BHEP}} test is performed with tuning parameter \eqn{\beta} chosen in dependence of the sample size n and the dimension d, namely \deqn{\beta=\frac{((2d+1)n/4)^(1/(d+4))}{\sqrt{2}}.}
#'
#' @references
#' Henze, N., Zirkler, B. (1990), A class of invariant consistent tests for multivariate normality, Commun.-Statist. - Th. Meth., 19:3595-3617, \href{https://doi.org/10.1080/03610929008830400}{DOI}
#'
#' @examples
#' test.HZ(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{HZ}}
#'
#' @export
test.HZ<-function(data, MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=HZ,tuning=NULL,repetitions = MC.rep)
  Testst=HZ(data)
  result <- list("Test" = "Henze-Zirkler", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Szekely-Rizzo (energy) test
#' @description
#' Performs the test of multivariate normality of Szekely and Rizzo (2005). Note that the scaled residuals use another scaling in the estimator of the covariance matrix!
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#' @param abb Stop criterium.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#'
#' @references
#' Szekely, G., Rizzo, M. (2005), A new test for multivariate normality, J. Multiv. Anal., 93:58-80, \href{https://doi.org/10.1016/j.jmva.2003.12.002}{DOI}
#'
#'
#' @examples
#' test.SR(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{SR}}
#'
#' @export
test.SR<-function(data, MC.rep=10000,alpha=0.05,abb=1e-8){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=SR,tuning=abb,repetitions = MC.rep)
  Testst=SR(data,abb)
  result <- list("Test" = "Szekely-Rizzo", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Manzotti-Quiroz test 1
#' @description
#' Performs the first test of multivariate normality of Manzotti and Quiroz (2001).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#'
#' @references
#' Manzotti, A., Quiroz, A.J. (2001), Spherical harmonics in quadratic forms for testing multivariate normality, Test, 10:87-104, \href{https://doi.org/10.1007/BF02595825}{DOI}
#'
#' @examples
#' test.MQ1(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=100)
#'
#' @seealso
#' \code{\link{MQ1}}
#'
#' @export
test.MQ1<-function(data, MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MQ1,tuning=NULL,repetitions = MC.rep)
  Testst=MQ1(data)
  result <- list("Test" = "Manzotti-Quiroz 1", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Manzotti-Quiroz test 2
#' @description
#' Performs the second test of multivariate normality of Manzotti and Quiroz (2001).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value
#' @param alpha level of significance of the test
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#'
#' @references
#' Manzotti, A., Quiroz, A.J. (2001), Spherical harmonics in quadratic forms for testing multivariate normality, Test, 10:87-104, \href{https://doi.org/10.1007/BF02595825}{DOI}
#'
#' @examples
#' test.MQ2(MASS::mvrnorm(50,c(0,1),diag(1,2)),MC.rep=500)
#'
#' @seealso
#' \code{\link{MQ2}}
#'
#' @export
test.MQ2<-function(data, MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=MQ2,tuning=NULL,repetitions = MC.rep)
  Testst=MQ2(data)
  result <- list("Test" = "Manzotti-Quiroz 2", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}


#' multivariate normality test of Cox and Small
#' @description
#' Performs the test of multivariate normality of Cox and Small (1978).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#' @param Points number of points to approximate the maximum functional on the unit sphere.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#'
#' @details
#' The test statistic is \eqn{T_{n,CS}=\max_{b\in\{x\in\mathbf{R}^d:\|x\|=1\}}\eta_n^2(b)},
#' where \deqn{\eta_n^2(b)=\frac{\left\|n^{-1}\sum_{j=1}^nY_{n,j}(b^\top Y_{n,j})^2\right\|^2-\left(n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^3\right)^2}{n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^4-1-\left(n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^3\right)^2}}.
#' Here, \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error. Note that the maximum functional has to be approximated by a discrete version, for details see Ebner (2012).
#'
#' @references
#' Cox, D.R., Small, N.J.H. (1978), Testing multivariate normality, Biometrika, 65:263-272.
#'
#' Ebner, B. (2012), Asymptotic theory for the test for multivariate normality by Cox and Small, Journal of Multivariate Analysis, 111:368-379.
#'
#' @examples
#' \donttest{test.CS(MASS::mvrnorm(10,c(0,1),diag(1,2)),MC.rep=100)}
#'
#' @seealso
#' \code{\link{CS}}
#'
#' @export
test.CS<-function(data, MC.rep=1000,alpha=0.05,Points=NULL){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=CS,tuning=Points,repetitions = MC.rep)
  Testst=CS(data,Points)
  result <- list("Test" = "Cox and Small", "param" = NULL, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Henze-Jimenes-Gamero test of multivariate normality
#'
#' Computes the multivariate normality test of Henze and Jimenes-Gamero (2019) in dependence of a tuning parameter \code{a}.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}#' }
#'
#' @references
#' Henze, N., Jimenez-Gamero, M.D. (2019) "A new class of tests for multinormality with i.i.d. and garch data based on the empirical moment generating function", TEST, 28, 499-521, \href{https://doi.org/10.1007/s11749-018-0589-z}{DOI}
#'
#' @examples
#' test.HJG(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=1.5,MC.rep=500)
#'
#' @seealso
#' \code{\link{HJG}}
#'
#' @export
test.HJG<-function(data,a=1,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=HJG,tuning=a,repetitions = MC.rep)
  Testst=HJG(data,a)
  result <- list("Test" = "Henze-Jimenes-Gamero", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}


#' The Henze-Visagie test of multivariate normality
#'
#' Computes the multivariate normality test of Henze and Visagie (2019).
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' Note that \code{a=Inf} returns the limiting test statistic with value \code{2*\link{MSkew} + \link{MRSSkew}}.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#' @references
#' Henze, N., Visagie, J. (2019) "Testing for normality in any dimension based on a partial differential equation involving the moment generating function", to appear in Ann. Inst. Stat. Math., \href{https://doi.org/10.1007/s10463-019-00720-8}{DOI}
#'
#' @examples
#' test.HV(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=5,MC.rep=500)
#' test.HV(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=Inf,MC.rep=500)
#'
#' @seealso
#' \code{\link{HV}}
#'
#' @export
test.HV<-function(data,a=5,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=HV,tuning=a,repetitions = MC.rep)
  Testst=HV(data,a)
  result <- list("Test" = "Henze-Visagie", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Henze-Jimenes-Gamero-Meintanis test of multivariate normality
#'
#' Computes the test statistic of the Henze-Jimenes-Gamero-Meintanis test.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#' @references
#' Henze, N., Jimenes-Gamero, M.D., Meintanis, S.G. (2019), Characterizations of multinormality and corresponding tests of fit, including for GARCH models, Econometric Th., 35:510-546, \href{https://doi.org/10.1017/S0266466618000154}{DOI}.
#'
#' @examples
#' \donttest{test.HJM(MASS::mvrnorm(10,c(0,1),diag(1,2)),a=2.5,MC=100)}
#'
#' @seealso
#' \code{\link{HJM}}
#'
#' @export
test.HJM<-function(data,a=1.5,MC.rep=500,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=HJM,tuning=a,repetitions = MC.rep)
  Testst=HJM(data,a)
  result <- list("Test" = "Henze-Jimenes-Gamero-Meintanis", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Pudelko test of multivariate normality
#'
#' Computes the (approximated) Pudelko test of multivariate normality.
#'
#' This functions evaluates the test statistic with the given data and the specified parameter \code{r}. Since since one has to calculate the supremum of a function inside a d-dimensional Ball of radius \code{r}. In this implementation the \code{\link{optim}} function is used.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#' @param r a positive number (radius of Ball)
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#' @references
#' Pudelko, J. (2005), On a new affine invariant and consistent test for multivariate normality, Probab. Math. Statist., 25:43-54.
#'
#' @examples
#' test.PU(MASS::mvrnorm(20,c(0,1),diag(1,2)),r=2,MC=100)
#'
#' @seealso
#' \code{\link{PU}}
#'
#' @export
test.PU<-function(data,MC.rep=10000,alpha=0.05,r=2){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=PU,tuning=r,repetitions = MC.rep)
  Testst=PU(data,r)
  result <- list("Test" = "Pudelko", "param" = r, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Doerr-Ebner-Henze test of multivariate normality based on harmonic oscillator
#'
#' Computes the multivariate normality test of Doerr, Ebner and Henze (2019) based on zeros of the harmonic oscillator.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#' @references
#' Doerr, P., Ebner, B., Henze, N. (2019) "Testing multivariate normality by zeros of the harmonic oscillator in characteristic function spaces" \href{https://arxiv.org/abs/1909.12624}{arXiv:1909.12624}
#'
#' @examples
#' test.DEHT(MASS::mvrnorm(20,c(0,1),diag(1,2)),a=1,MC=500)
#'
#' @seealso
#' \code{\link{DEHT}}
#'
#' @export
test.DEHT<-function(data,a=1,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=DEHT,tuning=a,repetitions = MC.rep)
  Testst=DEHT(data,a)
  result <- list("Test" = "DEH based on harmonic oscillator", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}



#' Doerr-Ebner-Henze test of multivariate normality based on a double estimation in a PDE
#'
#' Computes the multivariate normality test of Doerr, Ebner and Henze (2019) based on a double estimation in a PDE.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#' @references
#' Doerr, P., Ebner, B., Henze, N. (2019) "Testing multivariate normality by zeros of the harmonic oscillator in characteristic function spaces" \href{https://arxiv.org/abs/1909.12624}{arXiv:1909.12624}
#'
#' @examples
#' test.DEHU(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=1,MC=500)
#'
#' @seealso
#' \code{\link{DEHU}}
#'
#' @export
test.DEHU<-function(data,a=0.5,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=DEHU,tuning=a,repetitions = MC.rep)
  Testst=DEHU(data,a)
  result <- list("Test" = "DEH based on a double estimation in PDE", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}


#' Ebner-Henze-Strieder test of multivariate normality based on Fourier methods in a multivariate Stein equation
#'
#' Computes the multivariate normality test of Ebner, Henze and Strieder (2020) based on Fourier methods in a multivariate Stein equation.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#' @param MC.rep number of repetitions for the Monte Carlo simulation of the critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$Test}}{name of the test.}
#'         \item{\code{$param}}{value tuning parameter.}
#'         \item{\code{$Test.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'}
#'
#' @references
#' Ebner, B., Henze, N., Strieder, D. (2020) "Testing normality in any dimension by Fourier methods in a multivariate Stein equation" \href{https://arxiv.org/abs/2007.02596}{arXiv:2007.02596}
#'
#' @examples
#' test.EHS(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=1,MC=500)
#'
#' @seealso
#' \code{\link{EHS}}
#'
#' @export
test.EHS<-function(data,a=0.5,MC.rep=10000,alpha=0.05){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
  }else{warning("Wrong dimensions of data!")}}
  cv<-cv.quan(samplesize=n,dimension=d,quantile=1-alpha,statistic=EHS,tuning=a,repetitions = MC.rep)
  Testst=EHS(data,a)
  result <- list("Test" = "EHS Test", "param" = a, "Test.value" = Testst, "cv" = cv, "Decision" = (Testst > cv) )
  attr(result, "class") <- "mnt"
  return(result)
}
