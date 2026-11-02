# R/G2Gcov-package.R

#' G2Gcov: Discrete-Time Grassia(II)-Geometric Survival Models with Covariates
#'
#' Fits Grassia(II)-Geometric (G2G) discrete-time survival and churn models
#' with time-invariant or time-varying covariates via maximum likelihood
#' estimation.
#'
#' @keywords package
#' @importFrom stats aggregate optim as.formula ave model.matrix model.frame
#' @importFrom stats model.response terms delete.response .getXlevels qnorm
#' @importFrom stats dgamma qgamma
#' @importFrom graphics curve abline mtext
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end