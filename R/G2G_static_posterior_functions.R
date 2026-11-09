#' Covariate Effect for One Subject with Time-Invariant Covariates
#'
#' Internal helper shared by \code{\link{G2G_static_survival}},
#' \code{\link{G2G_static_posterior_p}} and
#' \code{\link{G2G_static_cond_survival}}. Codes one set of covariate values using
#' the terms and factor levels stored in the fit and returns the constant
#' covariate effect \eqn{\exp(x'\beta)}.
#'
#' @param solution A fitted model of class \code{"G2Gcov"} returned by
#'   \code{\link{G2G_static_MLE}}.
#' @param X_vect Covariate values for one subject, in any of the forms accepted
#'   by \code{X_vect} in \code{\link{G2G_static_survival}}. Only one set of values
#'   is accepted; anything else is an error.
#' @return A single number, \eqn{\exp(x'\beta)}.
#' @keywords internal
#' @importFrom stats model.frame model.matrix na.fail
G2G_static_effect <- function(solution, X_vect) {
  beta <- solution$par[-c(1:2)]

  # Build the one-row covariate matrix X from the covariate values supplied.
  if (is.data.frame(X_vect)) {
    # Original variable names; factors and dummies are coded as in the fit.
    if (nrow(X_vect) != 1L)
      stop("'X_vect' must have exactly one row. For several covariate ",
           "profiles, call this function once per profile (e.g. in a loop).",
           call. = FALSE)
    # na.fail so a missing value is an error instead of silently dropping the
    # row (which would return just S(0) = 1).
    mf <- model.frame(solution$terms, X_vect, xlev = solution$xlevels,
                      na.action = na.fail)
    X  <- model.matrix(solution$terms, mf)[, -1, drop = FALSE]
  } else {
    # Plain numbers for ONE subject: one value per coefficient (a single number
    # when the model has one coefficient), as a vector or a one-row matrix.
    # Several profiles are not accepted; the caller loops instead.
    if (is.matrix(X_vect)) {
      if (nrow(X_vect) != 1L)
        stop("'X_vect' must have exactly one row. For several covariate ",
             "profiles, call this function once per profile (e.g. in a loop).",
             call. = FALSE)
      X_vect <- X_vect[1, ]   # keeps the column names, which are used below
    }
    if (!is.numeric(X_vect))
      stop("'X_vect' must be a one-row data frame or numeric values.",
           call. = FALSE)
    if (length(beta) == 1L) {
      if (length(X_vect) != 1L)
        stop("The model has a single covariate coefficient, so 'X_vect' must ",
             "be a single number. For several values, call this function once ",
             "per value (e.g. in a loop).", call. = FALSE)
    } else if (length(X_vect) != length(beta)) {
      stop("'X_vect' must be a single row with ", length(beta), " values (one ",
           "per coefficient: ", paste(names(beta), collapse = ", "), "). For ",
           "several covariate profiles, call this function once per profile ",
           "(e.g. in a loop).", call. = FALSE)
    }
    # Named values (names = model-matrix column names) are matched to the
    # coefficients in any order; unnamed values are taken in coefficient order
    if (!is.null(names(X_vect))) {
      if (!setequal(names(X_vect), names(beta)))
        stop("Names of 'X_vect' must be: ",
             paste(names(beta), collapse = ", "), call. = FALSE)
      X_vect <- X_vect[names(beta)]          # reorder to match beta
    }
    X <- matrix(X_vect, nrow = 1)
  }

  # Constant covariate effect exp(x' * beta); c() turns the 1 x 1 matrix into a
  # number.
  c(exp(X %*% beta))
}

#' Check a Number of Periods
#'
#' Internal helper that stops unless \code{x} is a single whole number at or
#' above \code{min}.
#'
#' @param x Value to check.
#' @param name Argument name used in the error message.
#' @param min Smallest allowed value.
#' @return \code{NULL}, invisibly, when the value is valid.
#' @keywords internal
G2G_check_periods <- function(x, name, min) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      x < min || x != round(x))
    stop("'", name, "' must be a single whole number >= ", min, ".",
         call. = FALSE)
  invisible(NULL)
}

#' Posterior Mean of the Per-Period Event Probability Given History
#' (Time-Invariant Covariates)
#'
#' Updates the population distribution of the baseline rate using the number of
#' periods a subject has survived, and returns the subject's posterior mean
#' per-period event probability, which includes the effect of the subject's own
#' (constant) covariate values.
#'
#' With time-invariant covariates the subject's covariate effect is the constant
#' \eqn{c = \exp(x'\beta)}. Surviving \code{hist} periods gives
#' \eqn{C(\mathrm{hist}) = \mathrm{hist} \cdot c} and the posterior
#' \eqn{\lambda \sim \mathrm{Gamma}(r, \alpha + \mathrm{hist} \cdot c)}. The
#' subject's per-period event probability is \eqn{p = 1 - e^{-\lambda c}}, and
#' its posterior mean is
#' \deqn{1 - \left(\frac{\alpha + \mathrm{hist}\, c}
#'   {\alpha + (\mathrm{hist} + 1)\, c}\right)^r.}
#' This is also the probability that the subject has the event in the next
#' period given survival through \code{hist} periods, that is, 1 minus the first
#' step of \code{\link{G2G_static_cond_survival}} with \code{future = 1}. For a
#' subject whose covariates are all zero (\eqn{c = 1}) it is the baseline
#' probability \eqn{p_0}.
#'
#' @inheritParams G2G_static_survival
#' @param hist Number of periods the subject has survived without an event, a
#'   whole number >= 0. \code{hist = 0} returns the prior mean.
#' @return A single number: the posterior mean of the subject's per-period event
#'   probability given survival through \code{hist} periods.
#' @seealso \code{\link{G2G_static_cond_survival}},
#'   \code{\link{G2G_varying_posterior_p0}} for the baseline probability with
#'   time-varying covariates.
#' @export
#' @examples
#' \donttest{
#'   if (requireNamespace("survival", quietly = TRUE)) {
#'     library(survival)
#'     fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
#'     # A 60-year-old with Karnofsky score 70 who has survived 30 days
#'     G2G_static_posterior_p(fit, data.frame(age = 60, karno = 70), hist = 30)
#'   }
#' }
G2G_static_posterior_p <- function(solution, X_vect, hist) {
  G2G_check_periods(hist, "hist", 0)
  
  r     <- unname(solution$par[1])
  alpha <- unname(solution$par[2])
  
  # Constant covariate effect exp(x' * beta).
  C_c <- G2G_static_effect(solution, X_vect)
  
  # 1 - (b / (b + c))^r with b = alpha + hist * c, written to stay accurate when
  # c / b is small.
  b <- alpha + hist * C_c
  -expm1(-r * log1p(C_c / b))
}

#' Conditional Survival Over Future Periods Given History
#' (Time-Invariant Covariates)
#'
#' For a subject who has survived \code{hist} periods, computes the probability
#' of surviving each of the next \code{future} periods:
#' \eqn{P(T > \mathrm{hist} + k \mid T > \mathrm{hist})} for
#' \eqn{k = 1, \dots, \mathrm{future}}.
#'
#' With the constant covariate effect \eqn{c = \exp(x'\beta)}, surviving
#' \code{hist} periods updates the baseline rate to
#' \eqn{\mathrm{Gamma}(r, \alpha + \mathrm{hist} \cdot c)}, and averaging the
#' future survival over that posterior gives
#' \deqn{P(T > \mathrm{hist} + k \mid T > \mathrm{hist}) =
#'   \left(\frac{\alpha + \mathrm{hist}\, c}
#'   {\alpha + (\mathrm{hist} + k)\, c}\right)^r,}
#' which equals \eqn{S(\mathrm{hist} + k) / S(\mathrm{hist})} from
#' \code{\link{G2G_static_survival}}.
#'
#' @inheritParams G2G_static_survival
#' @param hist Number of periods the subject has survived without an event, a
#'   whole number >= 0. \code{hist = 0} gives the unconditional survival curve.
#' @param future Number of periods to forecast ahead, a whole number >= 1.
#' @return Numeric vector of length \code{future + 1}: a leading 1 (no additional
#'   periods), then the probability of surviving 1, 2, ..., \code{future}
#'   additional periods given survival through \code{hist} periods.
#' @seealso \code{\link{G2G_static_posterior_p}},
#'   \code{\link{G2G_varying_cond_survival}} for time-varying covariates.
#' @export
#' @examples
#' \donttest{
#'   if (requireNamespace("survival", quietly = TRUE)) {
#'     library(survival)
#'     fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
#'     # Survived 30 days; probability of surviving each of the next 10 days
#'     G2G_static_cond_survival(fit, data.frame(age = 60, karno = 70),
#'                              hist = 30, future = 10)
#'   }
#' }
G2G_static_cond_survival <- function(solution, X_vect, hist, future) {
  G2G_check_periods(hist,   "hist",   0)
  G2G_check_periods(future, "future", 1)

  r     <- unname(solution$par[1])
  alpha <- unname(solution$par[2])

  # Constant covariate effect exp(x' * beta).
  C_c <- G2G_static_effect(solution, X_vect)

  # Posterior rate b = alpha + hist * C; Delta_k = k * C is the added cumulative
  # effect over the next k periods.
  b     <- alpha + hist * C_c
  Delta <- C_c * seq_len(future)

  # (b / (b + Delta_k))^r, written to stay accurate when Delta_k / b is small.
  c(1, exp(-r * log1p(Delta / b)))
}
