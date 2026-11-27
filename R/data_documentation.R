# R/data_documentation.R

#' Kiwi Bubbles Dataset
#'
#' Person-period discrete-time trial-purchase dataset for "Kiwi Bubbles" (a
#' masked name for a shelf-stable multipack juice drink aimed primarily at
#' children) from Market 2 of an IRI BehaviorScan test market, observing 1,499
#' panelists over 24 weeks alongside weekly market-level marketing-mix
#' covariates.
#'
#' @format A data frame with 34,506 rows and 7 variables:
#' \describe{
#'   \item{id}{Integer panelist identifier.}
#'   \item{week}{Integer week index for the observation (running 1, 2, ...).}
#'   \item{num_weeks}{Total weeks observed for the panelist (integer).}
#'   \item{censor}{Censoring indicator on the panelist's final row: 1 if the
#'     panelist was censored without a trial purchase, 0 if a trial purchase
#'     occurred. 0 on all preceding rows.}
#'   \item{coupon}{Exponentially smoothed coupon stock variable capturing
#'     carryover effects of market-level coupon activity.}
#'   \item{anyp}{Market-level promotional activity measured as \%ACV with any
#'     promotion (the standard scanner-data sales-weighted promotional measure).}
#'   \item{status}{Discrete-time trial-purchase event indicator: 1 in the
#'     terminal period if the panelist made a first purchase, 0 otherwise.}
#' }
#'
#' @source IRI BehaviorScan "Kiwi Bubbles" test-market dataset (Market 2, first
#'   24 weeks), originally distributed at
#'   \url{https://brucehardie.com/datasets/kiwibubbles.zip} and processed via
#'   \file{data-raw/prepare_kb_data.R}.
#'
#' @examples
#' data(kb_data)
#' head(kb_data)
"kb_data"

#' Kiwi Bubbles Sample Dataset
#'
#' A 20-panelist subset of \code{\link{kb_data}} provided for fast unit
#' testing, minimal package checks, and vignette examples.
#'
#' @format A data frame with the same 7 variables as \code{\link{kb_data}}:
#' \describe{
#'   \item{id}{Integer panelist identifier.}
#'   \item{week}{Integer week index for the observation.}
#'   \item{num_weeks}{Total weeks observed for the panelist (integer).}
#'   \item{censor}{Censoring indicator (1 if censored on final row, 0 otherwise).}
#'   \item{coupon}{Exponentially smoothed coupon stock variable.}
#'   \item{anyp}{Market-level promotional activity (\%ACV with any promotion).}
#'   \item{status}{Discrete-time trial-purchase event indicator (1 if event
#'     occurred in the period, 0 otherwise).}
#' }
#'
#' @source Derived as a 20-panelist subset of \code{\link{kb_data}} via
#'   \file{data-raw/prepare_kb_data.R}.
#'
#' @examples
#' data(kb_data_sample)
#' head(kb_data_sample)
"kb_data_sample"
