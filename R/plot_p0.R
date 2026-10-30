#' Plot the baseline churn probability without covariates from a fitted G2Gcov Model
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
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{G2G_static_MLE}}, \code{\link{G2G_varying_MLE}}
#' @examples
#' \donttest{
#'   data(kb_data)
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
  
  col_dens <- "darkblue"; col_mean <- "red3"; col_med <- "darkorange"
  
  curve(dp0(x), from = 0.001, to = 0.999, n = 1000, lwd = 2, col = col_dens,
        main = "Distribution of Churn Probability without Covariates (p0)",
        xlab = "Churn Probability without Covariates (p0)", ylab = "Density")
  abline(v = mean_p0,   col = col_mean, lwd = 2, lty = 2)
  abline(v = median_p0, col = col_med,  lwd = 2, lty = 3)
  
  legend("topright",
         legend = c("Density",
                    sprintf("Mean = %.3f",   mean_p0),
                    sprintf("Median = %.3f", median_p0)),
         col = c(col_dens, col_mean, col_med),
         lwd = 2, lty = c(1, 2, 3),
         bty = "o", bg = "white", inset = 0.02)
  
  invisible(list(mean = mean_p0, median = median_p0))
}
