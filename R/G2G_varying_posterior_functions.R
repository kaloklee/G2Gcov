#' Per-Period Covariate Effects for a Covariate Path
#'
#' Internal helper shared by \code{\link{G2G_varying_survival}},
#' \code{\link{G2G_varying_posterior_p0}} and
#' \code{\link{G2G_varying_cond_survival}}. Codes a covariate path using the
#' terms and factor levels stored in the fit and returns the per-period
#' covariate effect \eqn{\exp(x_j'\beta)}, one value per row of the path.
#'
#' @param solution A fitted model of class \code{"G2Gcov"} returned by
#'   \code{\link{G2G_varying_MLE}}.
#' @param X_mat Covariate path with one row per period: a data frame holding the
#'   original covariate variables (factors and dummies are coded as in the fit;
#'   missing values are an error), or a matrix whose columns are named like the
#'   model coefficients (any order) or, if unnamed, already in coefficient order.
#'   \code{NULL} or a path with zero rows means "no periods" and returns
#'   \code{numeric(0)}.
#' @param arg Name of the argument being checked, used in error messages.
#' @return Numeric vector with one effect per period (\code{numeric(0)} for an
#'   empty path).
#' @keywords internal
#' @importFrom stats model.frame model.matrix na.fail
G2G_varying_effects <- function(solution, X_mat, arg = "X_mat") {
  beta <- solution$par[-(1:2)]

  # An empty path (no periods) contributes nothing.
  if (is.null(X_mat) || NROW(X_mat) == 0L)
    return(numeric(0))

  if (is.data.frame(X_mat)) {
    # Original variable names; factors and dummies coded as in the fit.
    # na.fail so a missing value cannot silently drop a period.
    mf <- model.frame(solution$terms, X_mat,
                      xlev = solution$xlevels, na.action = na.fail)
    X  <- model.matrix(solution$terms, mf)[, -1, drop = FALSE]
  } else {
    # Matrix input: use as given, after checking that it lines up with the
    # coefficients.
    X <- as.matrix(X_mat)
    # Named columns must match the coefficient names (any order, reordered below).
    if (!is.null(colnames(X))) {
      if (!setequal(colnames(X), names(beta)))
        stop("Column names of '", arg, "' must be: ",
             paste(names(beta), collapse = ", "), call. = FALSE)
      X <- X[, names(beta), drop = FALSE]      # reorder to match beta
    } else if (ncol(X) != length(beta)) {
      # Unnamed columns: only the number of columns can be checked, so they must
      # already be in the coefficient order.
      stop("'", arg, "' must have ", length(beta), " columns.", call. = FALSE)
    }
    if (anyNA(X))
      stop("'", arg, "' contains missing values.", call. = FALSE)
  }

  # Per-period covariate effect exp(x_t' * beta); c() turns the one-column
  # matrix into a vector.
  c(exp(X %*% beta))
}

#' Posterior Mean of the Baseline Per-Period Event Probability Given History
#'
#' Updates the population distribution of the baseline rate using a customer's
#' observed history and returns the posterior mean of the baseline per-period
#' event probability \eqn{p_0 = 1 - e^{-\lambda}}, i.e. the per-period event
#' probability at covariate value zero (the reference levels of any factors).
#'
#' A customer who survived \code{t} periods with covariate path
#' \eqn{z(1), \dots, z(t)} has posterior \eqn{\lambda \sim
#' \mathrm{Gamma}(r, \alpha + C(t))}, where \eqn{C(t) = \sum_{j=1}^{t}
#' e^{z(j)\beta'}}. The posterior mean of \eqn{p_0} is then
#' \deqn{1 - \left(\frac{\alpha + C(t)}{\alpha + C(t) + 1}\right)^r.}
#' The history enters only through \eqn{C(t)}.
#'
#' @inheritParams G2G_varying_survival
#' @param X_hist Covariate values for each period the customer has survived,
#'   one row per period (t = 1, ..., T), in any of the forms accepted by
#'   \code{X_mat} in \code{\link{G2G_varying_survival}}. The history must contain
#'   only event-free periods. \code{NULL} (no history) returns the prior mean.
#' @return A single number: the posterior mean of the baseline per-period event
#'   probability given survival through the full history.
#' @export
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE(
#'     "Surv(week, status) ~ coupon + anyp",
#'     data    = kb_data,
#'     subject = "id"
#'   )
#'   # A customer who has survived the first 8 weeks
#'   hist <- kb_data[kb_data$id == kb_data$id[1], c("coupon", "anyp")][1:8, ]
#'   G2G_varying_posterior_p0(fit, hist)
#' }
#'
G2G_varying_posterior_p0 <- function(solution, X_hist = NULL) {
  r     <- unname(solution$par[1])
  alpha <- unname(solution$par[2])

  # C(t): cumulative covariate effect over the history (0 if there is none).
  C_hist <- sum(G2G_varying_effects(solution, X_hist, "X_hist"))

  # 1 - (b / (b + 1))^r with b = alpha + C(t), written to stay accurate when
  # 1 / b is small.
  b <- alpha + C_hist
  -expm1(-r * log1p(1 / b))
}

#' Conditional Survival Over Planned Future Periods Given History
#'
#' For a customer who has survived the periods in a covariate history, computes
#' the probability of surviving each of the next \code{s} periods under a planned
#' covariate path: \eqn{P(T > t + k \mid T > t)} for \eqn{k = 1, \dots, s}.
#'
#' The customer's history updates the distribution of the baseline rate to
#' \eqn{\mathrm{Gamma}(r, \alpha + C(t))}, and averaging the future survival over
#' that posterior gives
#' \deqn{P(T > t + k \mid T > t) =
#'   \left(\frac{\alpha + C(t)}{\alpha + C(t + k)}\right)^r,}
#' where \eqn{C(m) = \sum_{j=1}^{m} e^{z(j)\beta'}} accumulates the history and
#' then the planned covariates. This equals \eqn{S(t + k) / S(t)} from
#' \code{\link{G2G_varying_survival}} for the combined path. The history length
#' \eqn{t} and the horizon \eqn{s} are taken from the number of rows of
#' \code{X_hist} and \code{X_future}.
#'
#' @inheritParams G2G_varying_survival
#' @param X_hist Covariate values for each period the customer has survived,
#'   one row per period, in any of the forms accepted by \code{X_mat} in
#'   \code{\link{G2G_varying_survival}}. The history must contain only
#'   event-free periods. \code{NULL} (no history) gives the unconditional
#'   survival curve for the planned path.
#' @param X_future Planned covariate values for the periods to forecast, one row
#'   per future period (periods t + 1, ..., t + s), in the same form as
#'   \code{X_hist}. At least one row is required.
#' @return Numeric vector of length s + 1: a leading 1 (no additional periods),
#'   then the probability of surviving 1, 2, ..., s additional periods given
#'   survival through the history.
#' @export
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE(
#'     "Surv(week, status) ~ coupon + anyp",
#'     data    = kb_data,
#'     subject = "id"
#'   )
#'   # Survived 8 weeks; plan 4 more weeks with a coupon and no other promotion
#'   hist   <- kb_data[kb_data$id == kb_data$id[1], c("coupon", "anyp")][1:8, ]
#'   future <- data.frame(coupon = rep(1, 4), anyp = rep(0, 4))
#'   G2G_varying_cond_survival(fit, hist, future)
#' }
#'
G2G_varying_cond_survival <- function(solution, X_hist, X_future) {
  r     <- unname(solution$par[1])
  alpha <- unname(solution$par[2])

  e_future <- G2G_varying_effects(solution, X_future, "X_future")
  if (length(e_future) == 0L)
    stop("'X_future' must have at least one row.", call. = FALSE)

  # C(t): cumulative effect over the history (0 if there is none).
  C_hist <- sum(G2G_varying_effects(solution, X_hist, "X_hist"))

  # Posterior rate b = alpha + C(t); Delta_k = C(t + k) - C(t) is the running
  # total of the planned effects.
  b     <- alpha + C_hist
  Delta <- cumsum(e_future)

  # (b / (b + Delta_k))^r, written to stay accurate when Delta_k / b is small.
  c(1, exp(-r * log1p(Delta / b)))
}
