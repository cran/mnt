#' Mardias measure of multivariate sample skewness
#' @description
#' This function computes the classical invariant measure of multivariate sample skewness due to Mardia (1970).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return value of sample skewness in the sense of Mardia.
#'
#' @details
#' Multivariate sample skewness due to Mardia (1970) is defined by
#' \deqn{b_{n,d}^{(1)}=\frac{1}{n^2}\sum_{j,k=1}^n(Y_{n,j}^\top Y_{n,k})^3,}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error. Note that for \eqn{d=1}, we have a measure proportional to the squared sample skewness.
#'
#' @references
#' Mardia, K.V. (1970), Measures of multivariate skewness and kurtosis with applications, Biometrika, 57:519–530.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467–506.
#'
#' @examples
#' MSkew(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
MSkew <-function(data)
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    S=((n-1)/n)*stats::var(data)
    Y=(data-mean(data))/sqrt(S)
    return((mean(Y^3))^2)
    }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    return(sum(Djk^3)/n^2)
  }
}



#' Mardias measure of multivariate sample kurtosis
#' @description
#' This function computes the classical invariant measure of multivariate sample kurtosis due to Mardia (1970).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return value of sample kurtosis in the sense of Mardia.
#'
#' @details
#' Multivariate sample kurtosis due to Mardia (1970) is defined by
#' \deqn{b_{n,d}^{(2)}=\frac{1}{n}\sum_{j=1}^n\|Y_{n,j}\|^4,}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}.To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' @references
#' Mardia, K.V. (1970), Measures of multivariate skewness and kurtosis with applications, Biometrika, 57:519–530.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467–506.
#'
#' @examples
#' MKurt(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
MKurt <-function(data)
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    S=((n-1)/n)*stats::var(data)
    Y=(data-mean(data))/sqrt(S)
    return(mean(Y^4))
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
  dmean=colMeans(data)
  normY=rep(0,n)
  Sn=((n-1)/n)*stats::cov(data)
  S=solve(Sn)
  for (i in 1:n)
  {
    normY[i]=(data[i,]-dmean)%*%S%*%(data[i,]-dmean)
  }
  return(mean(normY^2))
  }
}



#' Koziols measure of multivariate sample kurtosis
#' @description
#' This function computes the invariant measure of multivariate sample kurtosis due to Koziol (1989).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return value of sample kurtosis in the sense of Koziol.
#'
#' @details
#' Multivariate sample kurtosis due to Koziol (1989) is defined by
#' \deqn{\widetilde{b}_{n,d}^{(2)}=\frac{1}{n^2}\sum_{j,k=1}^n(Y_{n,j}^\top Y_{n,k})^4,}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error. Note that for \eqn{d=1}, we have a measure proportional to the squared sample kurtosis.
#'
#' @references
#' Koziol, J.A. (1989), A note on measures of multivariate kurtosis, Biom. J., 31:619–624.
#'
#' @examples
#' KKurt(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
KKurt <-function(data)
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    S=((n-1)/n)*stats::var(data)
    Y=(data-mean(data))/sqrt(S)
    return((mean(Y^4))^2)
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    return(mean(Djk^4))
  }
}





#' multivariate skewness in the sense of Malkovich and Afifi
#' @description
#' This function computes the invariant measure of multivariate sample skewness due to Malkovich and Afifi (1973).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param Points points for approximation of the maximum on the sphere. \code{Points=NULL} generates 1000 uniformly distributed Points on the d dimensional unit sphere.
#'
#' @return value of sample skewness in the sense of Malkovich and Afifi.
#'
#' @details
#' Multivariate sample skewness due to Malkovich and Afifi (1973) is defined by
#' \deqn{b_{n,d,M}^{(1)}=\max_{u\in \{x\in\mathbf{R}^d:\|x\|=1\}}\frac{\left(\frac{1}{n}\sum_{j=1}^n(u^\top X_j-u^\top \overline{X}_n )^3\right)^2}{(u^\top S_n u)^3},}
#' where \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' @references
#' Malkovich, J.F., and Afifi, A.A. (1973), On tests for multivariate normality, J. Amer. Statist. Ass., 68:176–179.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467–506.
#'
#' @examples
#' MASkew(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
MASkew<-function(data,Points = NULL)                           #Punkte sind hier Punkte auf der Einheitssph?re um das Maximum zu approximieren
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
    data=t(t(data))
    if(is.null(Points)){
      Unif= stats::runif(1000)
      Points = (Unif<=1/2)-(Unif>1/2)
    }
    GK=length(Points)
    Points=t(t(Points))
    bnm=rep(0,GK)
    mu=colMeans(data)
    Sn=((n-1)/n)*stats::cov(data)  #Ist hier die Normierung richtig im vgl. zur Funktion standard? (13.03.20) jetzt ja (16.03.20)
    for (j in 1:GK)
    {
      summe=0
      for (i in 1:n)
      {
        summe=summe+(t(Points[j,])%*%(data[i,]-mu))^3
      }
      summe=1/n*summe
      bnm[j]=summe^2/(t(Points[j,])%*%Sn%*%Points[j,])^3
    }
    return(max(bnm))
    }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
  if(is.null(Points)){
    Points = GVP(1000,d)
  }
  GK=dim(Points)[1]
  bnm=rep(0,GK)
  mu=colMeans(data)
  Sn=((n-1)/n)*stats::cov(data)  #Ist hier die Normierung richtig im vgl. zur Funktion standard? (13.03.20) jetzt ja (16.03.20)
  for (j in 1:GK)
  {
    summe=0
    for (i in 1:n)
    {
      summe=summe+(t(Points[j,])%*%(data[i,]-mu))^3
    }
    summe=1/n*summe
    bnm[j]=summe^2/(t(Points[j,])%*%Sn%*%Points[j,])^3
  }
  return(max(bnm))
  }
}




#' multivariate kurtosis in the sense of Malkovich and Afifi
#' @description
#' This function computes the invariant measure of multivariate sample kurtosis due to Malkovich and Afifi (1973).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param Points points for approximation of the maximum on the sphere. \code{Points=NULL} generates 1000 uniformly distributed Points on the d dimensional unit sphere.
#'
#' @return value of sample kurtosis in the sense of Malkovich and Afifi.
#'
#' @details
#' Multivariate sample skewness due to Malkovich and Afifi (1973) is defined by
#' \deqn{b_{n,d,M}^{(1)}=\max_{u\in \{x\in\mathbf{R}^d:\|x\|=1\}}\frac{\left(\frac{1}{n}\sum_{j=1}^n(u^\top X_j-u^\top \overline{X}_n )^3\right)^2}{(u^\top S_n u)^3},}
#' where \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' @references
#' Malkovich, J.F., and Afifi, A.A. (1973), On tests for multivariate normality, J. Amer. Statist. Ass., 68:176–179.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467–506.
#'
#' @examples
#' MAKurt(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
MAKurt<-function(data,Points = NULL)                           #Punkte sind hier Punkte auf der Einheitssph?re um das Maximum zu approximieren
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
    data=t(t(data))
    if(is.null(Points)){
      Unif= stats::runif(1000)
      Points = (Unif<=1/2)-(Unif>1/2)
    }
    GK=length(Points)
    Points=t(t(Points))
    bnm=rep(0,GK)
    mu=colMeans(data)
    Sn=((n-1)/n)*stats::cov(data)  #Ist hier die Normierung richtig im vgl. zur Funktion standard? (13.03.20) jetzt ja (16.03.20)
    for (j in 1:GK)
    {
      summe=0
      for (i in 1:n)
      {
        summe=summe+(t(Points[j,])%*%(data[i,]-mu))^4
      }
      summe=1/n*summe
      bnm[j]=summe/(t(Points[j,])%*%Sn%*%Points[j,])^2
    }
    return(max(bnm))
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
    if(is.null(Points)){
      Points = GVP(1000,d)
    }
    GK=dim(Points)[1]
    bnm=rep(0,GK)
    mu=colMeans(data)
    Sn=((n-1)/n)*stats::cov(data)  #Ist hier die Normierung richtig im vgl. zur Funktion standard? (13.03.20) jetzt ja (16.03.20)
    for (j in 1:GK)
    {
      summe=0
      for (i in 1:n)
      {
        summe=summe+(t(Points[j,])%*%(data[i,]-mu))^4
      }
      summe=1/n*summe
      bnm[j]=summe/(t(Points[j,])%*%Sn%*%Points[j,])^2
    }
    return(max(bnm))
  }
}



#' multivariate skewness of Móri, Rohatgi and Székely
#' @description
#' This function computes the invariant measure of multivariate sample skewness due to Móri, Rohatgi and Székely (1993).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return value of sample skewness in the sense of Móri, Rohatgi and Székely.
#'
#' @details
#' Multivariate sample skewness due to Móri, Rohatgi and Székely (1993) is defined by
#' \deqn{\widetilde{b}_{n,d}^{(1)}=\frac{1}{n}\sum_{j=1}^n\|Y_{n,j}\|^2\|Y_{n,k}\|^2Y_{n,j}^\top Y_{n,k},}
#' where \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error. Note that for \eqn{d=1}, it is equivalent to skewness in the sense of Mardia.
#'
#' @references
#' Móri, T. F., Rohatgi, V. K., Székely, G. J. (1993), On multivariate skewness and kurtosis, Theory of Probability and its Applications, 38:547–551.
#'
#' Henze, N. (2002), Invariant tests for multivariate normality: a critical review, Statistical Papers, 43:467–506.
#'
#' @export
MRSSkew<-function(data)
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){return(MSkew(data))}else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Yjk=matrix(0,n,n)
    return(mean(Rquad%*%t(Rquad)*Djk))
  }
}



#' Statistic of the BHEP-test
#' @description
#' This function returns the value of the statistic of the Baringhaus-Henze-Epps-Pulley (BHEP) test as in Henze and Wagner (1997).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param a positive numeric number (tuning parameter).
#'
#' @return value of the test statistic.
#'
#' @details
#' The test statistic is \deqn{BHEP_{n,\beta}=\frac{1}{n} \sum_{j,k=1}^n \exp\left(-\frac{\beta^2\|Y_{n,j}-Y_{n,k}\|^2}{2}\right)- \frac{2}{(1+\beta^2)^{d/2}} \sum_{j=1}^n \exp\left(- \frac{\beta^2\|Y_{n,j}\|^2}{2(1+\beta^2)} \right) + \frac{n}{(1+2\beta^2)^{d/2}}.}
#' Here, \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' @references
#' Henze, N., and Wagner, T. (1997), A new approach to the class of BHEP tests for multivariate normality, J. Multiv. Anal., 62:1–23, \href{https://doi.org/10.1006/jmva.1997.1684}{DOI}
#'
#' Epps T.W., Pulley L.B. (1983), A test for normality based on the empirical characteristic function, Biometrika, 70:723-726, \href{https://doi.org/10.1093/biomet/70.3.723}{DOI}
#'
#' @examples
#' BHEP(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
BHEP<-function(data,a=1){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    mu=mean(data)
    sn2=((n-1)/n)*stats::var(data)
    data=(data-mu)/sqrt(sn2)
    SUMME1=0
    SUMME2=0
    for (j in 1:n)
    {
      SUMME2=SUMME2+exp(-a^2*data[j]^2/(2*(1+a^2)))
      for (k in 1:n)
      {
        SUMME1=SUMME1+exp(-a^2*(data[j]-data[k])^2/2)
      }
    }
    ret=1/n*SUMME1-(2/sqrt(1+a^2))*SUMME2+n/sqrt(1+2*a^2)
    return(ret)
    }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else if(!a>0){warning("tuning parameter a>0 needed!")
  } else{
  data=standard(data)
  Djk=data%*%t(data)
  Rquad=diag(Djk)
  Dj=matrix(Rquad,n,n)
  Y=exp((-a^2/2)*(Dj-2*Djk+t(Dj)))
  Y2=exp((-a^2/(2*(1+a^2)))*Rquad)
  ret=sum(Y)/n-2*((1+a^2)^(-d/2))*sum(Y2)+((1+2*a^2)^(-d/2))*n
#  ret=sum(Y)/n**2-2*((1+a^2)^(-d/2))/n*sum(Y2)+((1+2*a^2)^(-d/2))  #(korrigiert 13.03.20)
  return(ret)
  }
}




#' statistic of the Székely-Rizzo test
#' @description
#' This function returns the value of the statistic of the test of multivariate normality (also called \emph{energy test}) as in Székely and Rizzo (2005). Note that the scaled residuals use another scaling in the estimator of the covariance matrix as the other functions of the package \code{mnt}!
#' It is equivalent to the function \code{\link[energy:mvnorm.test]{mvnorm.e}}.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param abb Stop criterium.
#'
#' @return value of the test statistic.
#'
#' @references
#' Székely, G., and Rizzo, M. (2005), A new test for multivariate normality, J. Multiv. Anal., 93:58–80, \href{https://doi.org/10.1016/j.jmva.2003.12.002}{DOI}
#'
#' @examples
#' SR(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#'@seealso
#'\code{\link[energy:mvnorm.test]{mvnorm.e}}
#'
#' @export
SR<-function(data,abb=1e-8){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    data=t(t(data))
    return(SR(data,abb))
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
  data=sqrt((n-1)/n)*standard(data)
  h=gamma((d+1)/2)                                                #Konstante
  g=2*h/gamma(d/2)                                                #Konstante
  RY=R(data)                                                      #Norm der skalierten Residuen
  S1=rep(0,n)                                                     #Vektor zur berechnung der ersten Summe
  RM=matrix(0,n,1)                                                #Matrix zur berechnung Der Taylorsumme bis zum Grad k
  LOG=rep(TRUE,n)
  LOG2=rep(FALSE,n)                                               #Logischer n Vektor
  RM[,1]=(RY)^(2)/(2) *gamma(3/2)/(gamma(1+ d/2))                 #Taylorsumme f?r k=0
  k=1                                                             #Grad k der Taylorsumme
  while(any(LOG!=LOG2)){                                          #Abbrechen fals alle LOG=TRUE
    Temp=rep(0,n)                                               #n Vektor
    Temp[LOG]=1/(factorial(k)*(-2)^(k))*(RY[LOG])^(2+2*k)/((2*k+1)*(2*k+2))*gamma(k+ 3/2)/(gamma(k+1+ d/2)) #Berechnung der Taylorsumme falls LOG=TRUE
    Temp[is.nan(Temp)]=0
    RM=c(RM,Temp)                                               #Hinzuf?gen von Temp zu RM
    dim(RM)=c(n,k+1)                                            #Dimensions Wiederherstellung
    LOG=abs(RM[,k+1])>abb                                       #falls Abbruchkriterium erreicht ab?ndern von LOG auf FALSE
    k=k+1                                                       #Grad erh?hen
  }
  Tay=rowSums(RM)                                                 #Taylorsummen Vektor
  if(any((Tay<0)==TRUE)) print(length(Tay[Tay<0]))
  Tay[Tay<0]=RM[(Tay<0),1]
  S1=sqrt(2)*h/gamma(d/2) + sqrt(2/pi)*h*Tay                      #Vektor zur berechnung der Ersten Summe
  Djk=data%*%t(data)                                              #Matrix der Machalanobis Abst?nde und Winkel
  Rquad=RY^2                                                      #Machalanobis Abst?nde bzw. R^2
  Dj=matrix(Rquad,n,n)                                            #Matrix von Rquad
  S2=Dj-2*Djk+t(Dj)                                               #Matrix zur Berechnung der zweiten Summe
  diag(S2)=rep(0,n)                                               #Diagonale auf 0 setzen da wegen Rundung negative werte auftreten
  Sum2=sum(sqrt(S2))                                              #komponentenweise Wurzel
  ret=2*sum(S1) - n*g - Sum2/n                                    #Teststatistik
  return(ret)
  }
}



#' Statistic of the Henze-Zirkler test
#' @description
#' This function returns the value of the statistic of the \code{\link{BHEP}} test as in Henze and Zirkler (1990). The difference to the \code{\link{BHEP}} test is in the choice of the tuning parameter \eqn{\beta}.
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return value of the test statistic.
#'
#' @details
#' A \code{\link{BHEP}} test is performed with tuning parameter \eqn{\beta} chosen in dependence of the sample size n and the dimension d, namely \deqn{\beta=\frac{((2d+1)n/4)^(1/(d+4))}{\sqrt{2}}.}
#'
#' @references
#' Henze, N., and Zirkler, B. (1990), A class of invariant consistent tests for multivariate normality, Commun.-Statist. – Th. Meth., 19:3595–3617, \href{https://doi.org/10.1080/03610929008830400}{DOI}
#'
#' @examples
#' HZ(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @seealso
#' \code{\link{BHEP}}
#'
#' @export
HZ<-function(data){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
    a=(((2*d+1)*n/4)^(1/(d+4)))/sqrt(2)
    return(BHEP(data,a))
    }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
  a=(((2*d+1)*n/4)^(1/(d+4)))/sqrt(2)
  return(BHEP(data,a))
  }
}




#' first statistic of Manzotti and Quiroz
#' @description
#' This function returns the value of the first statistic of Manzotti and Quiroz (2001).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return Value of the test statistic
#'
#' @references
#' Manzotti, A., and Quiroz, A.J. (2001), Spherical harmonics in quadratic forms for testing multivariate normality, Test, 10:87–104, \href{https://doi.org/10.1007/BF02595825}{DOI}
#'
#' @examples
#' MQ1(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
MQ1<-function(data){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){warning("Wrong dimensions of data!")} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else {
  data=standard(data)
  x=U(data)					#Winkelkomponente
  a2<-A2(d)
  a4<-A4(d)
  a5<-A5(d)
  a12<-A12(d)
  a13<-A13(d)
  h3<-H3(d)
  h6<-H6(d)
  h7<-H7(d)
  h9<-H9(d)
  h10<-H10(d)
  h15<-H15(d)
  m=choose(d+3,4)+choose(d+2,3)-1	#Anzahl der Funktionen
  vmat=matrix(0,m,n)			#m x n Matrix
  c=sqrt(norm(d,1))				#Konstante zur Normierung der Funktionen vom Grad 1
  for (j in 1:n){				#Spaltenweise Berechnung von vmat
    vmat[,j]=c(x[j,]*c,
               f2(x[j,],d,a2),
               f3(x[j,],d,h3),
               f4(x[j,],d,a4),
               f5(x[j,],d,a5,h3),
               f6(x[j,],d,a2,h6),
               f7(x[j,],d,h7),
               f8(x[j,],d,h7),		#harmonische Kugelfunktionen vom Grad 1 bis 4
               f9( x[j,],d,a2,h9),
               f10(x[j,],d,a4,h10),
               f11(x[j,],d,a5,h10,h3),
               f12(x[j,],d,a12),
               f13(x[j,],d,a13,h3),
               f14(x[j,],d,a4,h6),
               f15(x[j,],d,a5,h15))
  }
  v=rowMeans(vmat)*sqrt(n)		#Zeilen-Mittelwert von vmat und Normierung
  ret=t(v)%*%v				#quadratische Form
  return(as.numeric(ret))
  }
}



#' second statistic of Manzotti und Quiroz
#' @description
#' This function returns the value of the second statistic of Manzotti und Quiroz (2001).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#'
#' @return Value of the test statistic
#'
#' @references
#' Manzotti, A., and Quiroz, A.J. (2001), Spherical harmonics in quadratic forms for testing multivariate normality, Test, 10:87–104, \href{https://doi.org/10.1007/BF02595825}{DOI}
#'
#' @examples
#' MQ2(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
MQ2<-function(data){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){warning("Wrong dimensions of data!")} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else {
  data=standard(data)
  x=U(data)							#Winkelkomponente
  rad=R(data)							#radiale Komponente
  rad3=rad^3							#radiale Komponente hoch 3
  a2<-A2(d)
  a4<-A4(d)
  a5<-A5(d)
  a12<-A12(d)
  a13<-A13(d)
  h3<-H3(d)
  h6<-H6(d)
  h7<-H7(d)
  h9<-H9(d)
  h10<-H10(d)
  h15<-H15(d)				#Neuberechnung der globalen Hilfsmatrizen falls n?tig
  m=choose(d+1,2)+d+1					#Anzahl der Funktionen
  vmat=matrix(0,m,n)					#m x n Matrix
  V=matrix(0,m,m)						#Kovarianzmatrix V
  a=d*(d+2)*(d+4)						#Konstante a Erwarungswert von R^6
  b=gamma((d+1)/2)/gamma(d/2)				#Konstante b Erwarungswert von R
  diag(V)=a
  V[m-1,m-1]=d-2*b^2
  V[m,m-1]=V[m-1,m]=d*(d+2)-2*(d+1)*b^2		#V Definition
  V[m,m]=a-2*((d+1)^2)*(b^2)
  invV=solve(V)						#Inverse von V
  c=sqrt(norm(d,1))						#Konstante zur Normierung der Funktionen vom Grad 1
  for (j in 1:n){
    vmat[1:(m-2),j]=c(x[j,]*c,
                      f2(x[j,],d,a2),
                      f3(x[j,],d,h3))*rad3[j]	# harmonische Kugelfunktionen vom Grad 1 bis 2 mal R^3
    vmat[m-1,j]=rad[j]-sqrt(2)*b			# R-E(R)
    vmat[m,j]=rad3[j]-sqrt(2)*(d+1)*b		# R^3-E(R^3)
  }
  v=rowMeans(vmat)*sqrt(n)				#Zeilen-Mittelwert von vmat und Normierung
  ret=t(v)%*%invV%*%v					#quadratische Form
  return(as.numeric(ret))
  }
}



#' Statistic of the test of Cox and Small
#' @description
#' This function returns the (approximated) value of the test statistic of the test of Cox and Small (1978).
#'
#' @param data a n x d matrix of d dimensional data vectors.
#' @param Points points for approximation of the maximum on the sphere. \code{Points=NULL} generates 5000 uniformly distributed Points on the d dimensional unit sphere.
#'
#' @return approximation of the value of the test statistic of the test of Cox and Small (1978).
#'
#' @details
#' The test statistic is \eqn{T_{n,CS}=\max_{b\in\{x\in\mathbf{R}^d:\|x\|=1\}}\eta_n^2(b)},
#' where \deqn{\eta_n^2(b)=\frac{\left\|n^{-1}\sum_{j=1}^nY_{n,j}(b^\top Y_{n,j})^2\right\|^2-\left(n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^3\right)^2}{n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^4-1-\left(n^{-1}\sum_{j=1}^n(b^\top Y_{n,j})^3\right)^2}}.
#' Here, \eqn{Y_{n,j}=S_n^{-1/2}(X_j-\overline{X}_n)}, \eqn{j=1,\ldots,n}, are the scaled residuals, \eqn{\overline{X}_n} is the sample mean and \eqn{S_n} is the sample covariance matrix of the random vectors \eqn{X_1,\ldots,X_n}. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error. Note that the maximum functional has to be approximated by a discrete version, for details see Ebner (2012).
#'
#' @references
#' Cox, D.R. and Small, N.J.H. (1978), Testing multivariate normality, Biometrika, 65:263–272.
#'
#' Ebner, B. (2012), Asymptotic theory for the test for multivariate normality by Cox and Small, Journal of Multivariate Analysis, 111:368–379.
#'
#' @examples
#' CS(MASS::mvrnorm(50,c(0,1),diag(1,2)))
#'
#' @export
CS<-function(data,Points = NULL)
{
  n    = dim(data)[1]
  d    = dim(data)[2]
  if (is.null(n)){warning("Wrong dimensions of data!")} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else {
    if(is.null(Points)){
      Points = GVP(5000,d)
    }
  data=standard(data)
  GK=dim(Points)[1]
  etan=rep(0,GK)
  for (j in 1:GK)
  {
    Hilfsvar=H1(data,Points[j,])
    etan[j]=(Hilfsvar[1]-Hilfsvar[2])/(Hilfsvar[3]-1-Hilfsvar[2])
  }
  return(max(etan))
  }
}





#' Henze-Jiménes-Gamero test statistic
#'
#' Computes the test statistic of the Henze-Jimenes-Gamero test.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' @param data a n x d numeric matrix of data values.
#' @param a positive numeric number (tuning parameter).
#'
#' @return The value of the test statistic.
#'
#' @references
#' Henze, N., Jiménez-Gamero, M.D. (2019) "A new class of tests for multinormality with i.i.d. and garch data based on the empirical moment generating function", TEST, 28, 499-521, \href{https://doi.org/10.1007/s11749-018-0589-z}{DOI}
#'
#' @examples
#' HJG(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=5)
#'
#' @export
HJG <- function(data,a=5){
  n    = dim(data)[1]
  d    = dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    data=t(t(data))
    return(HJG(data,a))
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else if(!a>1){warning("tuning parameter a>1 needed!")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Yplus.norm=Dj+2*Djk+t(Dj)

    T1 = (a)^(-d/2)*sum(exp(Yplus.norm/(4*a)))/n
    T2 = n/((a-1)^(d/2))
    T3 = -2/((a-0.5)^(d/2))*sum(exp(Rquad/(4*a-2)))

    Tstat = T1+T2+T3
    return(Tstat)
  }
}



#' statistic of the Henze-Visagie test
#'
#' Computes the test statistic of the Henze-Visagie test.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' Note that \code{a=Inf} returns the limiting test statistic with value \code{2*\link{MSkew} + \link{MRSSkew}}.
#'
#' @param data a n x d numeric matrix of data values.
#' @param a numeric number greater than 1 (tuning parameter).
#'
#' @return The value of the test statistic.
#'
#' @references
#' Henze, N., Visagie, J. (2019) "Testing for normality in any dimension based on a partial differential equation involving the moment generating function", to appear in Ann. Inst. Stat. Math., \href{https://doi.org/10.1007/s10463-019-00720-8}{DOI}
#'
#' @examples
#' HV(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=5)
#' HV(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=Inf)
#'
#' @export
HV <- function(data,a=5){
  n    = dim(data)[1]
  d    = dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    data=t(t(data))
    return(HV(data,a))
    }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else if(!a>1){warning("tuning parameter a>1 needed!")
  } else if(a==Inf){
      return(2*MSkew(data) + MRSSkew(data))
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Yplus.norm=Dj+2*Djk+t(Dj)
    T1  = exp(Yplus.norm/(4*a))
    T2a = Djk
    T2b = Yplus.norm/(2*a)
    T2c = d/(2*a)
    T2d = Yplus.norm/(4*a^2)
    T2  = T2a-T2b+T2c+T2d
    T12 = sum(T1*T2)

    Tstat = (pi/a)^(d/2)*T12/n
    return(16*(a^(2+d/2))*Tstat/(pi)^(d/2))
  }
}



#' statistic of the Henze-Jiménes-Gamero-Meintanis test
#'
#' Computes the test statistic of the Henze-Jiménes-Gamero-Meintanis test.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the function returns an error.
#'
#' @param data a n x d numeric matrix of data values.
#' @param a positive numeric number (tuning parameter).
#'
#' @return The value of the test statistic.
#'
#' @references
#' Henze, N., Jiménes-Gamero, M.D., Meintanis, S.G. (2019), Characterizations of multinormality and corresponding tests of fit, including for GARCH models, Econometric Th., 35:510–546, \href{https://doi.org/10.1017/S0266466618000154}{DOI}.
#'
#' @examples
#' HJM(MASS::mvrnorm(20,c(0,1),diag(1,2)),a=2.5)
#'
#' @export
HJM <- function(data,a){
  n    = dim(data)[1]
  d    = dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    data=t(t(data))
    return(HJM(data,a))
    }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else if(!a>0){warning("tuning parameter a>0 needed!")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    help.exp=array(0,dim=c(n,n,n,n))
    for(j in 1:n)
    {
      dataj=matrix(data[j,],n,d,byrow=T)+data
      D1=dataj%*%t(dataj)
      for (l in 1:n)
      {
        datalp=matrix(data[l,],n,d,byrow=T)+data
        datalm=matrix(data[l,],n,d,byrow=T)-data
        D2p=datalp%*%t(datalp)
        D2m=datalm%*%t(datalm)
        Dkmp=dataj%*%t(datalp)
        Dkmm=dataj%*%t(datalm)
        Rquad1=diag(D1)
        Rquad2p=diag(D2p)
        Rquad2m=diag(D2m)
        Dj1=matrix(Rquad1,n,n)
        Dj2p=matrix(Rquad2p,n,n)
        Dj2m=matrix(Rquad2m,n,n)
        help.exp[j,l,,]=exp((Dj1-t(Dj2m))/(4*a))*cos(Dkmm/(2*a))+exp((Dj1-t(Dj2p))/(4*a))*cos(Dkmp/(2*a))
      }
    }
    help.sum1=sum(help.exp)
    help.sum2=sum(exp((Dj-t(Dj))/(4*a))*cos(Djk/(2*a)))
    Tstat = help.sum1/(2*n^3)-(2/n)*help.sum2+n
    return(Tstat)
  }
}



#' Statistic of the Pudelko test
#'
#' Approximates the test statistic of the Pudelko test.
#'
#' This functions evaluates the test statistic with the given data and the specified parameter \code{r}. Since since one has to calculate the supremum of a function inside a d-dimensional Ball of radius \code{r}. In this implementation the \code{\link{optim}} function is used.
#'
#' @param data a n x d numeric matrix of data values.
#' @param r a positive number (radius of Ball)
#'
#' @return approximate Value of the test statistic
#'
#' @references
#' Pudelko, J. (2005), On a new affine invariant and consistent test for multivariate normality, Probab. Math. Statist., 25:43–54.
#'
#' @examples
#' PU(MASS::mvrnorm(20,c(0,1),diag(1,2)),r=2)
#'
#' @export
PU<-function(data,r=2){
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){warning("Wrong dimensions of data!")} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else if(!r>0){warning("radius r>0 needed!")
  } else{
  data=standard(data)
  if(d>3) {it=1} else {it=5}
  vek=rep(0,it)
  if(it==1){start=U(t(MASS::mvrnorm(it,rep(0,d),diag(1,d))))*stats::runif(it)^2*r
  }else{
  start=U(MASS::mvrnorm(it,rep(0,d),diag(1,d)))*stats::runif(it)^2*r}
  for (i in 1:it) {
    vek[i]= stats::optim(start[i,],empchar,control=list(maxit=1000),data=data,r=r)$value
  }
  ret=-min(vek)*sqrt(n)
  return(ret)
}
}



#' Statistic of the DEH test based on harmonic oscillator
#'
#' Computes the test statistic of the DEH test.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a n x d numeric matrix of data values.
#' @param a positive numeric number (tuning parameter).
#'
#' @return The value of the test statistic.
#'
#' @references
#' Dörr, P., Ebner, B., Henze, N. (2019) "Testing multivariate normality by zeros of the harmonic oscillator in characteristic function spaces" \href{https://arxiv.org/abs/1909.12624}{arXiv:1909.12624}
#'
#' @examples
#' DEHT(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=1)
#'
#' @export
DEHT <- function(data,a=1) #DEHT
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
    mu=mean(data)
    sn2=((n-1)/n)*stats::var(data)
    data=(data-mu)/sqrt(sn2)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Y=Rquad%*%t(Rquad)*exp((-1/(4*a))*(Dj-2*Djk+t(Dj)))
    Y2=Rquad*(Rquad+2*d*a*(2*a+1))*exp((-1/(2*(1+2*a)))*Rquad)
    ret=(pi/a)^(d/2)*sum(Y)/n-((2*(2*pi)^(d/2))/((2*a+1)^(2+d/2)))*sum(Y2)+n*(pi^(d/2)/((a+1)^(2+d/2)))*(a*(a+1)*d^2+d*(d+2)/4)
    return(ret*d^(-2)*(a/pi)^(d/2))
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(a < 0 | !is.numeric(a)){
    warning("The parameter a needs to be a positive number! Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Y=Rquad%*%t(Rquad)*exp((-1/(4*a))*(Dj-2*Djk+t(Dj)))
    Y2=Rquad*(Rquad+2*d*a*(2*a+1))*exp((-1/(2*(1+2*a)))*Rquad)
    ret=(pi/a)^(d/2)*sum(Y)/n-((2*(2*pi)^(d/2))/((2*a+1)^(2+d/2)))*sum(Y2)+n*(pi^(d/2)/((a+1)^(2+d/2)))*(a*(a+1)*d^2+d*(d+2)/4)
    return(ret*d^(-2)*(a/pi)^(d/2))
  }
}



#' Statistic of the DEH test based on a double estimation in PDE
#'
#' Computes the test statistic of the DEH based on a double estimation in PDE test.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' @param data a (d,n) numeric matrix containing the data.
#' @param a positive numeric number (tuning parameter).
#'
#' @return The value of the test statistic.
#'
#' @references
#' Dörr, P., Ebner, B., Henze, N. (2019) "A new test of multivariate normality by a double estimation in a characterizing PDE" \href{https://arxiv.org/abs/1911.10955}{arXiv:1911.10955}
#'
#' @export
DEHU <- function(data,a)
{
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    n=length(data)
    d=1
    mu=mean(data)
    sn2=((n-1)/n)*stats::var(data)
    data=(data-mu)/sqrt(sn2)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Hilf=Dj-2*Djk+t(Dj)
    Y=(Rquad%*%t(Rquad) - (Dj+t(Dj))*(Hilf+2*d*a*(2*a-1))/(4*a^2)+(16*d^2*a^3*(a-1)+4*d*(d+2)*a^2+(8*d*a^2-(8+4*d)*a)*Hilf+Hilf^2)/(16*a^4))*exp(-Hilf/(4*a))
    ret=(pi/a)^(d/2)*sum(Y)/n
    return(ret*d^(-2)*(a/pi)^(d/2))
  }else{warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(a < 0 | !is.numeric(a)){
    warning("The parameter a needs to be a positive number! Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else{
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    Hilf=Dj-2*Djk+t(Dj)
    Y=(Rquad%*%t(Rquad) - (Dj+t(Dj))*(Hilf+2*d*a*(2*a-1))/(4*a^2)+(16*d^2*a^3*(a-1)+4*d*(d+2)*a^2+(8*d*a^2-(8+4*d)*a)*Hilf+Hilf^2)/(16*a^4))*exp(-Hilf/(4*a))
    ret=(pi/a)^(d/2)*sum(Y)/n
    return(ret*d^(-2)*(a/pi)^(d/2))
  }
}

#' Statistic of the EHS test based on a multivariate Stein equation
#'
#' Computes the test statistic of the EHS test based on a multivariate Stein equation.
#'
#' This functions evaluates the teststatistic with the given data and the specified tuning parameter \code{a}.
#' Each row of the data Matrix contains one of the n (multivariate) sample with dimension d. To ensure that the computation works properly
#' \eqn{n \ge d+1} is needed. If that is not the case the test returns an error.
#'
#' Note that \code{a=Inf} returns the limiting test statistic with value \code{2*\link{MSkew} + \link{MRSSkew}} and \code{a=0} returns the value of the limit statistic
#' \deqn{T_{n,0}=\frac{d}{2}-2^{\frac{d}{2}+1}\frac{1}{n}\sum_{j=1}^n\|Y_{n,j}\|^2\exp(-\frac{\|Y_{n,j}\|^2}{2}).}
#'
#' @param data a (d,n) numeric matrix containing the data.
#' @param a positive numeric number (tuning parameter).
#'
#' @return The value of the test statistic.
#'
#' @references
#' Ebner, B., Henze, N., Strieder, D. (2020) "Testing normality in any dimension by Fourier methods in a multivariate Stein equation" \href{https://arxiv.org/abs/2007.02596}{arXiv:2007.02596}
#'
#' @export
EHS<-function(data,a=1)
{
  require(mnt)
  n=dim(data)[1]
  d=dim(data)[2]
  if (is.null(n)){if (is.vector(data)){
    data=t(t(data))
    return(EHS(data,a))
  } else {warning("Wrong dimensions of data!")}} else if(n <= d){
    warning("Wrong data size! \n Maybe data needs to be transposed. Check the description for more Information.")
  } else if(!is.numeric(data)){
    stop("The data contains non-numeric entries! Check the description for more Information.")
  } else if(!a>=0){warning("tuning parameter a>0 needed!")
  } else if (a==0) {
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    return(d/2-2^(d/2+1)*mean(Rquad*exp(-Rquad/2)))
  } else if (a==Inf){return(2*MSkew(data)+MRSSkew(data))
  } else {
    data=standard(data)
    Djk=data%*%t(data)
    Rquad=diag(Djk)
    Dj=matrix(Rquad,n,n)
    SUM1=Djk*exp(-(1/(4*a))*(Dj-2*Djk+t(Dj)))
    SUM2=Rquad*exp(-Rquad/(4*a+2))/(2*a+1)
    ret=(pi/a)^(d/2)*sum(SUM1)/n-2*(2*pi/(2*a+1))^(d/2)*sum(SUM2)+n*(pi/(a+1))^(d/2)*d/(2*(a+1))
    return(ret)
  }
}
