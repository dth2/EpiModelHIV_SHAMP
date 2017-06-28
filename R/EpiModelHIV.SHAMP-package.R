#' HIV Transmission Dynamics among MSM and Heterosexuals
#'
#' \tabular{ll}{
#'    Package: \tab EpiModelHIV.SHAMP\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.0.0\cr
#'    Date: \tab 2017-06-27\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @details
#' EpiModelHIV.SHAMP provides extensions to EpiModel HIV to allow for simulating HIV transmission
#' dynamics among heterosexual and MSM of up to 5 race catagories simultaniously along with immigration.
#'
#' @name EpiModelHIV.SHAMP-package
#' @aliases EpiModelHIV.SHAMP
#'
#' @import EpiModel EpiModelHPC network networkDynamic tergmLite tergm ergm bindata
#' @importFrom stats rbinom rgeom rmultinom rpois runif simulate rnbinom plogis
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr group_by summarise
#'
#' @useDynLib EpiModelHIV.SHAMP
#'
#' @docType package
#' @keywords package msm het
#'
NULL
