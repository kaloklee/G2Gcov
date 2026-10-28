# R/data_documentation.R

#' KB dataset
#'
#' @description
#' Cleaned dataset created by the script `data-raw/prepare_kb_data.R`.
#' The raw input is read from `data-raw/kb_covars.csv` and column names
#' are converted to lower case for consistency.
#'
#' @format
#' A data frame with 34,506 rows (1,499 subjects observed for up to 24
#' weeks) and the following columns:
#' \describe{
#'   \item{id}{Integer subject identifier.}
#'   \item{week}{Integer week index for the observation.}
#'   \item{num_weeks}{Total weeks observed for the subject (integer).}
#'   \item{censor}{1 on the final row of a subject who had not yet had the
#'                  event by the end of observation; 0 on every other row.
#'                  A subject's event occurred in week \code{num_weeks}
#'                  unless that row has \code{censor = 1}.}
#'   \item{coupon}{Coupon stock covariate that measures the aggregate effects
#'                  of distributing coupons to potential customers.}
#'   \item{anyp}{Promotions covariate that measures the percent of stores in 
#'                the market with trade promotions within a given week.}
#' }
#'
#' @source
#' Internal raw CSV: \file{data-raw/kb_covars.csv}, processed by
#' \file{data-raw/prepare_kb_data.R}. 
#'
#' @examples
#' head(kb_data)
"kb_data"

#' KB Sample Dataset
#'
#' @description
#' A small subset of \code{\link{kb_data}} (all 50 rows of the first 22
#' subjects) provided for use in quick examples and tests where the full
#' dataset would be unnecessarily slow.
#'
#' @format
#' A data frame with the same columns as \code{kb_data}:
#' \describe{
#'   \item{id}{Integer subject identifier.}
#'   \item{week}{Integer week index for the observation.}
#'   \item{num_weeks}{Total weeks observed for the subject (integer).}
#'   \item{censor}{1 on the final row of a censored subject; 0 otherwise.}
#'   \item{coupon}{Coupon stock covariate.}
#'   \item{anyp}{Trade promotions covariate.}
#' }
#'
#' @source
#' Derived from \code{\link{kb_data}}; see \file{data-raw/prepare_kb_data.R}.
"kb_covars_sample"
