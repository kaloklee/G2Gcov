# R/data_documentation.R

#' KB Dataset
#'
#' Person-period discrete-time survival dataset observing 1,499 subjects for up
#' to 24 weeks with marketing covariates.
#'
#' @format A data frame with 34,506 rows and 7 variables:
#' \describe{
#'   \item{id}{Integer subject identifier.}
#'   \item{week}{Integer week index for the observation (running 1, 2, ...).}
#'   \item{num_weeks}{Total weeks observed for the subject (integer).}
#'   \item{censor}{Censoring indicator on the subject's final row: 1 if the
#'     subject was censored before experiencing the event, 0 if the event
#'     occurred. 0 on all preceding rows.}
#'   \item{coupon}{Coupon stock covariate measuring aggregate distributed
#'     coupons.}
#'   \item{anyp}{Trade promotions covariate measuring percentage of stores with
#'     active promotions in the market during the week.}
#'   \item{status}{Discrete-time failure event indicator: 1 in the terminal
#'     period if the subject experienced the event, 0 otherwise.}
#' }
#'
#' @source Simulated customer churn data modeled after retail cohort dynamics.
#'   Processed via \file{data-raw/prepare_kb_data.R}.
#'
#' @examples
#' data(kb_data)
#' head(kb_data)
"kb_data"

#' KB Sample Dataset
#'
#' A 20-subject random sample of \code{\link{kb_data}} provided for fast unit
#' testing, minimal package checks, and vignette examples.
#'
#' @format A data frame with the same 7 variables as \code{\link{kb_data}}:
#' \describe{
#'   \item{id}{Integer subject identifier.}
#'   \item{week}{Integer week index for the observation.}
#'   \item{num_weeks}{Total weeks observed for the subject (integer).}
#'   \item{censor}{Censoring indicator (1 if censored on final row, 0 otherwise).}
#'   \item{coupon}{Coupon stock covariate.}
#'   \item{anyp}{Trade promotions covariate.}
#'   \item{status}{Discrete-time failure event indicator (1 if event occurred in
#'     the period, 0 otherwise).}
#' }
#'
#' @source Derived as a random sample from \code{\link{kb_data}} via
#'   \file{data-raw/prepare_kb_data.R}.
#'
#' @examples
#' data(kb_data_sample)
#' head(kb_data_sample)
"kb_data_sample"