# MonteCarlo

#Quantiles

test.quantiles<-function(samplesize=c(20,50,100), dimension=1:3, quantiles=0.95, statistic, tuning=NULL, repetitions=100000)
{
  if(!is.null(tuning)){
  Erg=array(0,dim=c(repetitions,length(dimension),length(tuning),length(samplesize)))
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    for (l in 1:repetitions)
    {
    if (d==1) {X=stats::rnorm(max(samplesize))}else{X=MASS::mvrnorm(max(samplesize),rep(0,d),diag(1,d))}
    j.n=0
    for (n in samplesize)
    {
      j.n=j.n+1
      j.a=0
      for(a in tuning)
      {
        j.a=j.a+1
        if (d==1) {Erg[l,j.d,j.a,j.n]=statistic(X[1:n],a)}else{Erg[l,j.d,j.a,j.n]=statistic(X[1:n,],a)}
      }
    }
    }
  }
  help.quantiles=array(0,dim=c(length(dimension),length(tuning),length(samplesize), length(quantiles)))
  j.d=0
  for (d in dimension)
    {
      j.d=j.d+1
      j.n=0
      for (n in samplesize)
      {
        j.n=j.n+1
        j.a=0
        for(a in tuning)
        {
          j.a=j.a+1
          help.quantiles[j.d,j.a,j.n,]=sort(Erg[,j.d,j.a,j.n])[floor(quantiles*repetitions)]
        }
      }
  }
  }else{
    Erg=array(0,dim=c(repetitions,length(dimension),length(samplesize)))
    j.d=0
    for (d in dimension)
    {
      j.d=j.d+1
      for (l in 1:repetitions)
      {
        if (d==1) {X=stats::rnorm(max(samplesize))}else{X=MASS::mvrnorm(max(samplesize),rep(0,d),diag(1,d))}
        j.n=0
        for (n in samplesize)
        {
          j.n=j.n+1
          if (d==1) {Erg[l,j.d,j.n]=statistic(X[1:n])}else{Erg[l,j.d,j.n]=statistic(X[1:n,])}
        }
      }
    }
    help.quantiles=array(0,dim=c(length(dimension),length(samplesize), length(quantiles)))
    j.d=0
    for (d in dimension)
    {
      j.d=j.d+1
      j.n=0
      for (n in samplesize)
      {
        j.n=j.n+1
        help.quantiles[j.d,j.n,]=sort(Erg[,j.d,j.n])[floor(quantiles*repetitions)]
      }
   }
  }
  return(help.quantiles)
}


#' Monte Carlo simulation of quantiles for normality tests
#' @description
#' This function returns the quantiles of a test statistic with optional tuning parameter.
#'
#' @param samplesize samplesize for which the empirical quantile should be calculated.
#' @param dimension a natural number to specify the dimension of the multivariate normal distribution
#' @param quantile a number between 0 and 1 to specify the quantile of the empirical distribution of the considered test
#' @param statistic a function specifying the test statistic.
#' @param tuning the tuning parameter of the test statistic.
#' @param repetitions number of Monte Carlo runs.
#'
#' @return empirical quantile of the test statistic.
#'
#' @examples
#' cv.quan(samplesize=10, dimension=2,quantile=0.95, statistic=BHEP, tuning=2.5, repetitions=1000)
#'
#'@export
cv.quan<-function(samplesize, dimension, quantile, statistic, tuning=NULL, repetitions=100000)
{
  pb <- utils::txtProgressBar(style = 3)
  Erg=rep(0,repetitions)
  if(!is.null(tuning)){
    for (l in 1:repetitions)
    {
      utils::setTxtProgressBar(pb, l/repetitions)
      X=MASS::mvrnorm(samplesize,rep(0,dimension),diag(1,dimension))
      Erg[l]=statistic(X,tuning)
    }
  }else{
    for (l in 1:repetitions)
    {
      utils::setTxtProgressBar(pb, l/repetitions)
      X=MASS::mvrnorm(samplesize,rep(0,dimension),diag(1,dimension))
      Erg[l]=statistic(X)
    }
  }
  close(pb)
  return(sort(Erg)[floor(quantile*repetitions)])
}

#usethis::use_data(Ergebnisse.Quantile095, overwrite =T)
