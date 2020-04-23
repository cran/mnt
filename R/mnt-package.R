# ## usethis namespace: start
# #' @useDynLib mnt
# #' @importFrom Rcpp sourceCpp
# ## usethis namespace: end
# NULL


# Dim<<-0    #Setzen der Globalen variable Dimension eventuell benötigt für Manzotti


#' Print method for tests of multivariate normality
#' @description
#' Printing objects of class "mnt".
#'
#' @param x object of class "mnt".
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' A \code{mnt} object is a named list of numbers and character string, supplemented with \code{test} (the name of the teststatistic). \code{test} is displayed as a title.
#' The remaining elements are given in an aligned "name = value" format.
#'
#' @return
#' the argument x, invisibly, as for all \link{print} methods.
#'
#' @examples
#' print(test.DEHU(MASS::mvrnorm(50,c(0,1),diag(1,2)),a=1,MC=500))
#'
#' @export
print.mnt <- function(x, ...){
  cat("\n")
  cat("------------------------------------------------------------------------- \n")
  cat("\n")
  cat("         Test for multivariate normality with the", x$Test, " teststatistic.\n"  )
  cat("\n")
  if(!is.null(x$param)){cat("tuning parameter =", x$param," \n")}
  cat( x$Test, " = ", x$Test.value, " \n")
  cat("critical value =  ", x$cv, " (via monte carlo) \n")
  cat("\n")
  cat("\n")
  cat("------------------------------------------------------------------------- \n")
}



#' Simulated empirical 90\% quantiles of the tests contained in package \code{mnt}
#'
#' A dataset containing the empirical 0.9 quantiles of the tests for the dimensions \code{d=2,3,5} and samplesizes \code{n=20,50,100} based on a Monte Carlo Simulation study with 100000 repetitions. The following parameters were used:
#' \itemize{
#'   \item For \code{\link{BHEP}} the parameter \code{a=1},
#'   \item for \code{\link{HV}} the parameter \code{a=5},
#'   \item for \code{\link{HJG}} the parameter \code{a=1.5},
#'   \item for \code{\link{HJM}} the parameter \code{a=1.5},
#'   \item for \code{\link{DEHT}} the parameter \code{a=0.25},
#'   \item for \code{\link{DEHU}} the parameter \code{a=0.5},
#'   \item for \code{\link{CS}} the parameter \code{Points=NULL},
#'   \item for \code{\link{PU}} the parameter \code{r=2},
#'   \item for \code{\link{MASkew}} the parameter \code{Points=NULL},
#'   \item for \code{\link{MAKurt}} the parameter \code{Points=NULL},
#' }
#'
#' @format A data frame with 9 rows and 20 columns.
"Quantile09"



#' Simulated empirical 95\% quantiles of the tests contained in package \code{mnt}
#'
#' A dataset containing the empirical 0.95 quantiles of the tests for the dimensions \code{d=2,3,5} and samplesizes \code{n=20,50,100} based on a Monte Carlo Simulation study with 100000 repetitions. The following parameters were used:
#' \itemize{
#'   \item For \code{\link{BHEP}} the parameter \code{a=1},
#'   \item for \code{\link{HV}} the parameter \code{a=5},
#'   \item for \code{\link{HJG}} the parameter \code{a=1.5},
#'   \item for \code{\link{HJM}} the parameter \code{a=1.5},
#'   \item for \code{\link{DEHT}} the parameter \code{a=0.25},
#'   \item for \code{\link{DEHU}} the parameter \code{a=0.5},
#'   \item for \code{\link{CS}} the parameter \code{Points=NULL},
#'   \item for \code{\link{PU}} the parameter \code{r=2},
#'   \item for \code{\link{MASkew}} the parameter \code{Points=NULL},
#'   \item for \code{\link{MAKurt}} the parameter \code{Points=NULL},
#' }
#'
#' @format A data frame with 9 rows and 20 columns.
"Quantile095"


#' Simulated empirical 99\% quantiles of the tests contained in package \code{mnt}
#'
#' A dataset containing the empirical 0.99 quantiles of the tests for the dimensions \code{d=2,3,5} and samplesizes \code{n=20,50,100} based on a Monte Carlo Simulation study with 100000 repetitions. The following parameters were used:
#' \itemize{
#'   \item For \code{\link{BHEP}} the parameter \code{a=1},
#'   \item for \code{\link{HV}} the parameter \code{a=5},
#'   \item for \code{\link{HJG}} the parameter \code{a=1.5},
#'   \item for \code{\link{HJM}} the parameter \code{a=1.5},
#'   \item for \code{\link{DEHT}} the parameter \code{a=0.25},
#'   \item for \code{\link{DEHU}} the parameter \code{a=0.5},
#'   \item for \code{\link{CS}} the parameter \code{Points=NULL},
#'   \item for \code{\link{PU}} the parameter \code{r=2},
#'   \item for \code{\link{MASkew}} the parameter \code{Points=NULL},
#'   \item for \code{\link{MAKurt}} the parameter \code{Points=NULL},
#' }
#'
#' @format A data frame with 9 rows and 20 columns.
"Quantile099"
