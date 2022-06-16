#' @title Synthetics Data for Small Area Estimation using Hierarchical Bayesian Method under Zero Inflated Negative Binomial Distribution with non-sampled areas
#'
#' @description Dataset to simulate Small Area Estimation using Hierarchical Bayesian Method under Zero Inflated Negative Binomial Distribution with non-sampled areas
#'
#' This data contains \code{NA} values that indicates no sampled at one or more small areas. It uses the \code{\link{dataZINB}} with the direct estimates and the related variances in 5 small areas are missing.
#'
#'@usage data(dataZINBNS)
#'
#' @format A data frame with 50 rows and 4 variables::
#' \describe{
#'   \item{y}{Direct Estimation of y}
#'   \item{x1}{Auxiliary variable of x1}
#'   \item{x2}{Auxiliary variable of x2}
#'   \item{vardir}{Sampling Variance of y}
#' }
#'
"dataZINBNS"
