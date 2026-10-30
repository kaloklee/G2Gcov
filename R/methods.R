#' Print a Fitted G2Gcov Model
#'
#' @description
#' Displays a formatted summary of a G2G model fit returned by
#' \code{\link{G2G_static_MLE}} (time-invariant covariates) or
#' \code{\link{G2G_varying_MLE}} (time-varying covariates), including
#' convergence status, log-likelihood, AIC, and a parameter table with
#' standard errors and 95\% Wald confidence intervals.
#'
#' @details
#' Both fitting functions return an object of class \code{"G2Gcov"}, and this
#' one method prints either. It only relies on the elements the two have in
#' common: \code{par}, \code{par_stderr}, \code{par_lower}, \code{par_upper},
#' \code{value} and \code{convergence}. If any of the standard errors or
#' confidence limits are missing or \code{NA} (the fitting functions return
#' \code{NA} when the Hessian at the optimum is singular or not positive
#' definite), they print as \code{NA} and a note is added below the table.
#'
#' @param x A fitted model object of class \code{"G2Gcov"}.
#' @param digits Integer. Number of decimal places to print. Default 4.
#' @param ... Further arguments passed to or from other methods (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{G2G_static_MLE}}, \code{\link{G2G_varying_MLE}}
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
#'                           data = kb_data, subject = "id")
#'   print(fit)
#' }
#' @export
print.G2Gcov <- function(x, digits = 4, ...) {
  k <- length(x$par)

  # A fit from either G2G_static_MLE() or G2G_varying_MLE() is printed from the
  # same elements. Anything a fit does not carry is shown as NA instead of
  # making the print fail.
  or_na <- function(el) {
    v <- x[[el]]
    if (is.null(v)) rep(NA_real_, k) else v
  }
  par_names <- names(x$par)
  if (is.null(par_names)) par_names <- paste0("par", seq_len(k))
  value <- if (is.null(x$value)) NA_real_ else x$value

  conv_msg <- if (is.null(x$convergence)) {
    "Unknown"
  } else if (x$convergence == 0L) {
    "Converged (code 0)"
  } else {
    paste0("Did NOT converge (code ", x$convergence,
           ") -- interpret estimates with caution")
  }

  aic <- round(2 * k + 2 * value, digits)
  ll  <- round(-value, digits)

  cat("G2G Model Fit\n")
  cat(strrep("-", 52), "\n")
  cat(sprintf("%-18s %s\n", "Convergence:",    conv_msg))
  cat(sprintf("%-18s %s\n", "Log-likelihood:", ll))
  cat(sprintf("%-18s %s\n", "AIC:",            aic))
  cat(sprintf("%-18s %d\n", "Parameters:",     k))
  cat("\nCoefficients:\n")

  se  <- or_na("par_stderr")
  tab <- data.frame(
    Estimate  = round(x$par,           digits),
    Std.Error = round(se,              digits),
    Lower.95  = round(or_na("par_lower"), digits),
    Upper.95  = round(or_na("par_upper"), digits),
    row.names = par_names,
    check.names = FALSE
  )
  print(tab)

  if (anyNA(se))
    cat("\nNote: some standard errors are NA (the Hessian at the optimum is",
        "singular or not positive definite).\n")

  invisible(x)
}


#' Plot the Baseline Churn Probability Distribution of a Fitted G2Gcov Model
#'
#' @description
#' Plots the distribution of the baseline churn probability \eqn{p_0}, with all
#' covariates equal to zero, implied by a model fitted with
#' \code{\link{G2G_static_MLE}} (time-invariant covariates) or
#' \code{\link{G2G_varying_MLE}} (time-varying covariates). The mean and median
#' of \eqn{p_0} are marked with vertical lines and shown in a second title.
#'
#' @details
#' In the G2G model, \eqn{p_0 = 1 - e^{-\lambda}}, where \eqn{\lambda} follows a
#' Gamma distribution with shape \eqn{r} and rate \eqn{\alpha}. The fitted
#' \code{r} and \code{alpha} are on the original covariate scale, so the curve
#' describes a subject whose covariates are all zero, which can lie outside the
#' range of the data. The density is drawn exactly (by a change of variables,
#' not by simulation), and the mean and median have closed forms:
#' \eqn{1 - (\alpha / (1 + \alpha))^r} and
#' \eqn{1 - e^{-q}}, where \eqn{q} is the median of the Gamma distribution.
#' When \eqn{r < 1} the density is unbounded at 0, so the curve rises steeply
#' there and the mean lies above the median.
#'
#' @param solution A fitted model object of class \code{"G2Gcov"}.
#'
#' @return Invisibly returns a list with the \code{mean} and \code{median} of
#'   \eqn{p_0}.
#'
#' @seealso \code{\link{G2G_static_MLE}}, \code{\link{G2G_varying_MLE}}
#' @importFrom graphics curve abline mtext
#' @importFrom stats dgamma qgamma
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
#'                           data = kb_data, subject = "id")
#'   plot_p0(fit)
#' }
#' @export
plot_p0 <- function(solution) {
  r     <- unname(solution$par[1])
  alpha <- unname(solution$par[2])
  
  # Exact density of p0 = 1 - exp(-lambda), where lambda ~ Gamma(shape = r, rate = alpha)
  dp0 <- function(p) dgamma(-log1p(-p), shape = r, rate = alpha) / (1 - p)
  
  # Closed-form mean and median of p0
  mean_p0   <- 1 - (alpha / (1 + alpha))^r
  median_p0 <- 1 - exp(-qgamma(0.5, shape = r, rate = alpha))
  
  curve(dp0, from = 0.001, to = 0.999, n = 1000, lwd = 2, col = "darkblue",
        main = "Distribution of Baseline Churn Probability without Covariates (p0)",
        xlab = "Baseline Churn Probability (p0)", ylab = "Density")
  abline(v = mean_p0,   col = "red3",       lwd = 2, lty = 2)
  abline(v = median_p0, col = "darkorange", lwd = 2, lty = 3)
  
  # Second title: mean and median (line styles match the vertical lines)
  mtext(sprintf("Mean (dashed) = %.3f   |   Median (dotted) = %.3f", mean_p0, median_p0),
        side = 3, line = 0.4, cex = 0.9)
  
  invisible(list(mean = mean_p0, median = median_p0))
}

