#eventuell benoetigte Funktionen

#' L^2 random field
#'
#' Generates a L^2 random field
#'
#' This function calculates the critical values given a covariance kernel.
#'
#' @param d a positive integer indicating the dimension.
#' @param n.iter a positive integer indicating the number of iterations.
#' @param n.points a positive integer specifing ???
#' @param a Positive numeric number.
#' @param kov a with d^(-2) normed covariance kernel.
#'
#' @return A number (the critical value).
# @export
L2.random.field<-function(d, n.iter=1000, n.points=100, a, kov)          #Methode zur Berechnung der kritischen Werte bei gegebenem Kovarianzkern "kov" (Normierung mit d^(-2) beachten!)
{
  data = MASS::mvrnorm(n.points, rep(0,d), 1/(2*a) * diag(1,d))
  n = dim(data)[1]
  help = matrix(0, nrow=n, ncol=n)
  for (j in 1:n)
  {
    for (k in 1:n)
    {
      help[j,k] = kov( data[j,], data[k,])
    }
  }
  Z = MASS::mvrnorm(n.iter,rep(0,n),help)
  mean.erg = rep(0,n.iter)
  for (j in 1:n.iter)
  {
    mean.erg[j] = mean(Z[j,]^2)*d^(-2)                                    #Normierung!
  }
  return( sort( mean.erg)[ floor( 0.95 * n.iter)])
}


monte.carlo.p.value <- function(value, statistic ,d,n , n.iter = 10000, a){
  #    data = MASS::mvrnorm(n, rep(0,d), 1/(2*a) * diag(1,d))
  q = numeric(n.iter)
  pb <- utils::txtProgressBar(title="Example progress bar", label="0% done", min=0, max=100, initial=0,style = 3)
  for (i in 1:n.iter)
  {
    data = MASS::mvrnorm(n, rep(0,d),  diag(1,d))
    if(is.null(a)){ q[i] = statistic(data)
    } else{ q[i] = statistic(data,a)}
    utils::setTxtProgressBar(pb,i/n.iter*100)
  }
  return( 1-stats::ecdf(q)(value))
}

#' Helpfunction for covariance kernel
#'
#' A function that returns
#'
#' @param x a numeric vector.
#'
#' @return A number.
#'
# @examples
# hilf.e(1:20)
#'
#' @family helpfunctions
#'
# @export
hilf.e<-function(x)                                                #Hilfsfunktion Kovarianzkerne
{
  return( exp( -sum( x * x) / 2))
}

#' Euclidian norm
#'
#' A function that returns the euclidian norm of a vector.
#'
#' @param x a numeric vector.
#'
#' @return A number.
#'
# @examples
#  ENorm(1:20)
#'
#' @family helpfunctions
#'
# @export
ENorm<-function(x){
  sqrt( sum( x^2))
  }


#' Standarising data
#'
#' A function that return the standarized data.
#'
#' @param data a n x d matrix. ???
#'
#' @return A n x d matrix.
#'
# #' @example ...
#'
#' @family helpfunctions
#'
#' @export
standard<-function(data)                                                # Standardisierung
{
  n=dim(data)[1]
  d=dim(data)[2]
  Sn=((n-1)/n)*stats::cov(data)                                              # Empirische Kovarianzmatrix (richtige Normierung!)
  Snew = pracma::sqrtm(Sn)$Binv                                        # Berechnung von S_n^(-1/2),
  meanX=colMeans(data)                                                # Berechnung des empirischen Mittelwertes der Stichprobe
  Y=matrix(rep(0,d*n),ncol=d,nrow=n)                                  # Initialisierung der transformierten Daten Y_n,j
  for (i in 1:n)                                                      # Datentransformation
  {
    Y[i,]= Snew%*%(data[i,]-meanX)
  }
  return(Y)
}



# @export
#print.mnt <- function(x, ...){ #alt
#   cat("------------------------------------------------------------------------- \n")
#   cat("\n")
#   cat("         One-sample test for normality with the DEHT teststatistic.\n"  )
#   cat("\n")
#   cat("DEHT = ", x$Test.value, " \n")
#   cat("critical value =  ", x$cv, " \n")
#   cat("\n")
#   if (x$Decision == FALSE) {
#     cat("The test has no objection, that the sample is not  distributed. \n")
#   } else{
#     cat("The test rejects the assumption, that the sample is  distributed. \n")
#   }
#   cat("\n")
#   cat("------------------------------------------------------------------------- \n")
# }

#' Berechnet das negative des gewichteten Abstandes zwischen der empirischen und theoretischen charakteristischen Funktion an der Stelle x
#' @param	x d-dimensionaler Vektor
#' @param	data n x d Matrix der skalierten Residuen
#' @param r reelle Konstante des PUTests
#' @return reeller Wert
empchar=function(x,data,r){
  normx=ENorm(x) #Norm von x
  if(normx>r){ret=normx-r} #Ab?ndern der Funktion ausserhalb von r auf positive Werte
  else {n=dim(data)[1]
  y=t(t(x)%*%t(data)) #Vektor der Skalarprodukte
  ret=-abs((sum(exp(1i*y))/n-exp(-(normx^2)/2))/normx) # - Abstand (empirischen - theoretische char. Fnk. ), 1/normx Gewichtsfunktion
  }
  return(ret)
}

GVP<-function(Anzahl,d)   # Berechnung der Gleichverteilungsdiskretisierung
{
  mu=rep(0,d)
  E=diag(1,d,d)
  X2= MASS::mvrnorm(Anzahl,mu,E)
  GVPAO=matrix(rep(0,Anzahl*d),nrow=Anzahl,ncol=d)
  for (j in 1:Anzahl)
  {
    GVPAO[j,]=X2[j,]/(ENorm(X2[j,]))
  }
  return(GVPAO)
}

#' Hilfsfunktionen zur Berechnung der Teststatistik von Cox und Small
#' @param	data n x d Matrix der skalierten Residuen
#' @param vector vector on the unit sphere
#' @return reeller Wert
H1<-function(data,vector)
{
  n=dim(data)[1]
  d=dim(data)[2]
  summe1=rep(0,d)
  summe2=0
  summe3=0
  for (j in 1:n)
  {
    sp=sum(vector*data[j,])
    summe1=summe1+data[j,]*sp^2
    summe2=summe2+sp^3
    summe3=summe3+sp^4
  }
  summe1=1/n*summe1
  summe2=1/n*summe2
  summe3=1/n*summe3
  return(c((ENorm(summe1))^2,summe2^2,summe3))
}
