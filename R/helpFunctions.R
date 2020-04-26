#eventuell benoetigte Funktionen


ENorm<-function(x){
  sqrt( sum( x^2))
  }


#' Empirical scaled residuals
#'
#' A function that computes the scaled residuals of the data.
#'
#' @param data a n x d matrix of d dimensional data vectors..
#'
#' @return A n x d matrix of the scaled residuals.
#'
#' @export
standard<-function(data)                                                # Standardisierung
{
  n=dim(data)[1]
  d=dim(data)[2]
  Sn=((n-1)/n)*stats::cov(data)
  Snew = pracma::sqrtm(Sn)$Binv
  meanX=colMeans(data)
  Y=matrix(rep(0,d*n),ncol=d,nrow=n)
  for (i in 1:n)
  {
    Y[i,]= Snew%*%(data[i,]-meanX)
  }
  return(Y)
}


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
