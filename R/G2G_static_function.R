#' Fit G2G Model with Static Covariates
#'
#' @param par Parameter vector
#' @param df Data frame with time, status, and covariates
#' @return Negative log-likelihood value
#' @keywords internal
G2G_static_LL <- function(par, df) {
  r <- exp(par[1])
  alpha <- exp(par[2])
  beta <- par[-c(1:2)]
  
  y <- df[, 1]
  status <- df[, 2]
  X <- as.matrix(df[, -c(1:2)])
  
  # Uncensored likelihood
  uncen <- y[which(status == 1)]
  C_u <- exp(X[which(status == 1), ] %*% beta)
  LL_uncen <- logdiffexp(-r * log(1 + C_u * (uncen - 1) / alpha), 
                         -r * log(1 + C_u * uncen / alpha))
  
  # Censored likelihood
  cen <- y[which(status == 0)]
  C_c <- exp(X[which(status == 0), ] %*% beta)
  LL_cen <- -r * log(1 + C_c * cen / alpha)
  
  return(-sum(LL_uncen) - sum(LL_cen))
}

#' Fit a G2G Model with Time-Invariant Covariates
#'
#' @description
#' Fits a Grassia(II)-Geometric (G2G) survival model with time-invariant
#' covariates using maximum likelihood estimation via BFGS. Each subject is
#' represented by a single row. The survivor function at discrete period
#' \eqn{t} is
#' \eqn{S(t) = (1 + t \cdot \exp(\mathbf{x}'\boldsymbol{\beta}) /
#' \alpha)^{-r}},
#' where \eqn{r} and \eqn{\alpha} are the shape and rate of the Gamma
#' mixing distribution and \eqn{\boldsymbol{\beta}} are covariate
#' log-hazard coefficients.
#'
#' @param formula A formula of the form
#'   \code{Surv(time, status) ~ x1 + x2 + ...} where \code{time} is the
#'   observed period, \code{status} is the event indicator (1 = event,
#'   0 = censored), and \code{x1}, \code{x2}, etc. are time-invariant
#'   covariates.
#' @param data A data frame containing all variables referenced in
#'   \code{formula}. One row per subject.
#'
#' @return A list containing the optimization results from
#'   \code{\link[stats]{optim}} with the following additional elements:
#'   \item{par}{Named numeric vector of parameter estimates. The first two
#'     elements are \code{r (shape)} and \code{alpha (rate)}, returned on
#'     the natural (exponentiated) scale. Remaining elements are covariate
#'     log-hazard coefficients.}
#'   \item{par_stderr}{Standard errors derived from the observed Hessian,
#'     on the same scale as \code{par} (delta method for \code{r} and
#'     \code{alpha}). \code{NA} when the Hessian is singular or
#'     near-singular.}
#'   \item{par_upper}{Upper bounds of 95\% Wald confidence intervals. For
#'     \code{r} and \code{alpha} the interval is computed on the log scale
#'     and exponentiated, so it is not symmetric around the estimate.}
#'   \item{par_lower}{Lower bounds of 95\% Wald confidence intervals.}
#'   \item{value}{Negative log-likelihood at the optimum.}
#'   \item{convergence}{Convergence code from \code{optim} (0 = success).}
#'   \item{hessian}{Observed Hessian matrix at the optimum.}
#'
#' @seealso \code{\link{G2G_varying_MLE}} for the time-varying extension;
#'   \code{\link{plot_p0}} to visualize the fitted baseline propensity
#'   distribution.
#'
#' @importFrom stats optim model.matrix model.frame model.response
#' @examples
#' \donttest{
#'   if (requireNamespace("survival", quietly = TRUE)) {
#'     library(survival)
#'     fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
#'     print(fit)
#'   }
#' }
#' @export
G2G_static_MLE <- function(formula, data) {

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  vars <- all.vars(formula)
  missing_cols <- setdiff(vars, names(data))
  if (length(missing_cols) > 0L)
    stop("Column(s) not found in 'data': ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)

  na_cols <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1L))]
  if (length(na_cols) > 0L)
    stop("NA values found in column(s): ",
         paste(na_cols, collapse = ", "),
         ". Remove or impute missing values before fitting.", call. = FALSE)
  # ---------------------------------------------------------------------------

  # Read time and status straight from the columns named on the left-hand
  # side, so Surv() is never evaluated and the survival package need not be
  # attached.
  formula <- as.formula(formula)
  lhs     <- all.vars(formula[[2]])
  time    <- data[[lhs[1]]]
  status  <- data[[lhs[2]]]
  X <- model.matrix(stats::delete.response(stats::terms(formula)),
                    data = data)[, -1, drop = FALSE]

  if (!all(status %in% c(0, 1)))
    stop("Status column '", lhs[2], "' must contain only 0 and 1.",
         call. = FALSE)
  if (any(time < 1))
    stop("Time column '", lhs[1], "' must contain values >= 1.",
         call. = FALSE)
  if (sum(status) == 0L)
    stop("No events found (all status = 0). At least one event is required.",
         call. = FALSE)
  n_par <- 2L + ncol(X)
  if (nrow(data) <= n_par)
    stop("Insufficient data: need more rows than parameters (", n_par, ").",
         call. = FALSE)

  # Create data frame for optimization
  df <- data.frame(time = time, status = status, X)

  # Run optimization — muffle expected NaNs produced during BFGS exploration
  solution <- withCallingHandlers(
    optim(par     = c(log(0.5), log(0.5), rep(0, ncol(X))),
          fn      = G2G_static_LL,
          df      = df,
          method  = "BFGS",
          hessian = TRUE),
    warning = function(w) {
      if (grepl("NaNs produced", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )
  # Construct parameter names
  par_names <- c("r (shape)", "alpha (rate)", colnames(X))

  names(solution$par) <- par_names
  

  # Calculate standard errors and confidence intervals
  hess_inv <- tryCatch(solve(solution$hessian),
                       error = function(e) NULL)
  solution$par_stderr <- if (is.null(hess_inv)) {
    rep(NA_real_, length(solution$par))
  } else {
    sqrt(diag(hess_inv))
  }
  solution <- .natural_scale(solution)
  class(solution) <- c("G2Gcov", "list")
  return(solution)
}
