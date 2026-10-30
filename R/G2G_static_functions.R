#' Log-Likelihood for G2G Model with Static Covariates
#'
#' Internal function that computes the negative log-likelihood for the G2G model
#' with time-invariant covariates (one row per subject). Used by the
#' optimization routine in \code{\link{G2G_static_MLE}}.
#'
#' @param par Numeric vector of parameters on the optimizer's internal scale.
#'   The first two elements are log(r) and log(alpha) (both are exponentiated
#'   inside the function), and the remaining elements are the covariate
#'   coefficients.
#' @param df Data frame with columns, in this order: time (the observed period,
#'   a whole number >= 1), status (1 = event, 0 = censored), then the
#'   covariates. One row per subject.
#' @return Scalar value of the negative log-likelihood
#' @keywords internal
G2G_static_LL <- function(par, df) {
  # Overview: with time-invariant covariates each subject has one row and a constant
  # covariate effect C = exp(x' * beta). The survival function implied by the G2G
  # model is
  #   S(t) = (1 + t * C / alpha)^(-r),   with S(0) = 1.
  # Each subject contributes to the likelihood:
  #   - event in period t:   log[ S(t-1) - S(t) ]
  #   - censored at time t:  log[ S(t) ]
  # The function returns the NEGATIVE of the summed log-likelihood so that optim()
  # can minimize it.
  #
  # Expected layout of df: column 1 = time, column 2 = status (1 = event,
  # 0 = censored), columns 3 onward = covariates; one row per subject.
  # par = c(log r, log alpha, covariate coefficients).
  #
  # Back-transform r and alpha from the log scale used by the optimizer
  # (guarantees both are positive); the rest of par are the coefficients.
  r <- exp(par[1])
  alpha <- exp(par[2])
  beta <- par[-c(1:2)]
  
  # Unpack the data frame: observed period, event indicator and covariate matrix.
  y <- df[, 1]
  status <- df[, 2]
  X <- as.matrix(df[, -c(1:2)])
  
  # Uncensored likelihood
  # Event periods for subjects who had the event.
  uncen <- y[which(status == 1)]
  # Covariate effect exp(x' * beta) for those subjects. drop = FALSE keeps a
  # matrix when there is only one covariate (otherwise the subset collapses to a
  # vector and the matrix product fails with "non-conformable arguments").
  C_u <- exp(X[which(status == 1), , drop = FALSE] %*% beta)
  # log[ S(t-1) - S(t) ] for each subject, via logdiffexp(a, b) = log(exp(a) - exp(b))
  # computed stably on the log scale (defined in utils.R). For t = 1, S(0) = 1
  # falls out of the same formula.
  LL_uncen <- logdiffexp(-r * log(1 + C_u * (uncen - 1) / alpha), 
                         -r * log(1 + C_u * uncen / alpha))
  
  # Censored likelihood
  # Last observed period for censored subjects.
  cen <- y[which(status == 0)]
  # Covariate effect exp(x' * beta) for those subjects (drop = FALSE as above).
  C_c <- exp(X[which(status == 0), , drop = FALSE] %*% beta)
  # A censored subject survived through period t, so the contribution is
  # log S(t) = -r * log(1 + t * C / alpha).
  LL_cen <- -r * log(1 + C_c * cen / alpha)
  
  # Negative log-likelihood (optim() minimizes, so the sign is flipped).
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
#' Covariates are centered and scaled internally to keep the likelihood
#' numerically stable, and the results are converted back to the original
#' covariate scale.
#'
#' @param formula A formula of the form
#'   \code{Surv(time, status) ~ x1 + x2 + ...} where \code{time} is the
#'   observed period, \code{status} is the event indicator (1 = event,
#'   0 = censored), and \code{x1}, \code{x2}, etc. are time-invariant
#'   covariates. Transformations, interactions and factors are allowed on the
#'   right-hand side, but the intercept must be kept (no \code{- 1} or
#'   \code{+ 0}): the intercept column is always dropped, because \eqn{\alpha}
#'   plays its role.
#' @param data A data frame containing all variables referenced in
#'   \code{formula}. One row per subject. The observed times must be whole
#'   numbers >= 1 (discrete periods).
#'
#' @return A list of class \code{"G2Gcov"} containing the optimization results
#'   from \code{\link[stats]{optim}} with the following additional elements.
#'   All estimates, standard errors and confidence bounds are reported on the
#'   original covariate scale (only \code{hessian} is on the internal
#'   parameterization):
#'   \item{par}{Named numeric vector of parameter estimates. The first two
#'     elements are \code{r (shape)} and \code{alpha (rate)}, returned on
#'     the natural (exponentiated) scale. Remaining elements are covariate
#'     log-hazard coefficients in original units.}
#'   \item{par_stderr}{Standard errors on the original scale (named like
#'     \code{par}). For \code{r} and \code{alpha} these come from the delta
#'     method (e.g. SE(r) is approximately r * SE(log r)). If the Hessian at the
#'     optimum is not positive definite (or cannot be inverted) a warning is
#'     issued and these are \code{NA}.}
#'   \item{par_upper}{Upper bounds of 95\% Wald confidence intervals. For
#'     \code{r} and \code{alpha} the interval is formed on the log scale and
#'     exponentiated.}
#'   \item{par_lower}{Lower bounds of 95\% Wald confidence intervals.}
#'   \item{value}{Negative log-likelihood at the optimum.}
#'   \item{convergence}{Convergence code from \code{optim} (0 = success; a
#'     warning is issued for any other code).}
#'   \item{hessian}{Observed Hessian matrix at the optimum, on the internal
#'     parameterization (log r, log alpha, standardized coefficients), with
#'     row and column names.}
#'   \item{center}{Covariate means used internally for centering.}
#'   \item{scale}{Covariate standard deviations used internally for scaling.}
#'   \item{terms}{The \code{terms} object for the covariates (response
#'     removed), used to code new data in \code{\link{G2G_static_survival}}.}
#'   \item{xlevels}{Factor levels of the covariates at fitting time, also used
#'     to code new data.}
#'
#' @seealso \code{\link{G2G_varying_MLE}} for the time-varying extension;
#'   \code{\link{G2G_static_survival}} and \code{\link{G2G_static_retention}}
#'   for the fitted curves; \code{\link{plot_p0}} to visualize the fitted
#'   baseline propensity distribution.
#'
#' @importFrom stats optim model.matrix model.frame model.response terms
#' @importFrom stats delete.response
#' @importFrom stats .getXlevels
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

  # Workflow: (1) validate the inputs, (2) build and standardize the covariate
  # matrix, (3) maximize the likelihood with BFGS, (4) get standard errors and 95%
  # Wald intervals from the inverse Hessian, and (5) convert the results back to
  # the original covariate scale.
  #
  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)

  # Every variable named in the formula (time, status and covariates).
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

  # The formula is evaluated by model.frame()/model.matrix(), so transformations,
  # interactions and factor() terms are honored here.
  # Extract variables from formula
  mf <- model.frame(formula, data = data)
  # Surv response: column 1 = time, column 2 = status.
  y  <- model.response(mf)
  if (!inherits(y, "Surv") || !identical(attr(y, "type"), "right"))
    stop("The response must be a right-censored survival object, e.g. ",
         "Surv(time, status).", call. = FALSE)
  # The intercept column is dropped below, so it has to exist; without it the
  # first real covariate would be dropped instead.
  if (attr(attr(mf, "terms"), "intercept") == 0L)
    stop("The formula must keep its intercept (do not use '- 1' or '+ 0'): ",
         "alpha plays the role of the intercept.", call. = FALSE)
  # Design matrix without the intercept column (alpha plays the role of the
  # intercept); factors expand to dummy columns.
  X  <- model.matrix(formula, data = data)[, -1, drop = FALSE]

  # The model is for discrete periods 1, 2, 3, ... (the likelihood uses S(t - 1)),
  # so observed times must be whole numbers >= 1.
  if (any(y[, 1] < 1))
    stop("Observed times must be >= 1.", call. = FALSE)
  if (any(abs(y[, 1] - round(y[, 1])) > sqrt(.Machine$double.eps)))
    stop("Observed times must be whole numbers (discrete periods 1, 2, 3, ...).",
         call. = FALSE)

  # Sanity checks: at least one event, and more rows than parameters
  # (r, alpha and one coefficient per design-matrix column).
  if (sum(y[, 2]) == 0L)
    stop("No events found (all status = 0). At least one event is required.",
         call. = FALSE)
  n_par <- 2L + ncol(X)
  if (nrow(data) <= n_par)
    stop("Insufficient data: need more rows than parameters (", n_par, ").",
         call. = FALSE)

  # Standardize the covariates (center and scale) so exp(X %*% beta) cannot
  # overflow and the optimizer is well conditioned. This is an exact
  # reparameterization of the model: centering is absorbed by alpha and scaling
  # rescales the coefficients. Both are undone below.
  Xs  <- scale(X)
  ctr <- attr(Xs, "scaled:center")
  scl <- attr(Xs, "scaled:scale")
  if (any(scl == 0))
    stop("Covariate(s) with no variation: ",
         paste(names(scl)[scl == 0], collapse = ", "),
         ". Remove them before fitting.", call. = FALSE)

  # Create data frame for optimization
  # Layout expected by G2G_static_LL(): time, status, then the standardized covariates.
  df <- data.frame(time = y[, 1], status = y[, 2], Xs)

  # Run optimization — muffle expected NaNs produced during BFGS exploration
  # Starting values: r = alpha = 0.5 (passed to optim as log(0.5), because the
  # parameters are on the log scale) and all coefficients 0 (no covariate
  # effect). BFGS is unconstrained, which works because r and alpha are
  # optimized on the log scale. hessian = TRUE returns the Hessian of the negative
  # log-likelihood, used below for standard errors. The handler hides "NaNs
  # produced" warnings from extreme trial values; other warnings pass through.
  # At this point the parameters are on the internal scale: log r, log alpha and
  # coefficients for the standardized covariates.
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
  # optim() signals non-convergence (for example, reaching maxit) only through
  # this code, so surface it as a warning.
  if (solution$convergence != 0)
    warning("optim() did not converge (convergence code ", solution$convergence,
            "); estimates and standard errors may be unreliable.", call. = FALSE)

  # Construct parameter names
  par_names <- c("r (shape)", "alpha (rate)", colnames(Xs))

  names(solution$par) <- par_names
  # The Hessian stays on the internal parameterization, so its labels say so.
  hess_names <- c("log(r)", "log(alpha)", colnames(Xs))
  dimnames(solution$hessian) <- list(hess_names, hess_names)

  # Calculate standard errors and confidence intervals
  # Covariance matrix of the internal-scale estimates: the inverse of the Hessian
  # of the negative log-likelihood. The inverse is only used when the Hessian is
  # finite and positive definite (all eigenvalues > 0), i.e. the solution is a
  # proper minimum. Otherwise the standard errors are returned as NA with a
  # warning.
  hess_ok <- all(is.finite(solution$hessian)) &&
    min(eigen(solution$hessian, symmetric = TRUE, only.values = TRUE)$values) > 0
  hess_inv <- if (hess_ok) {
    tryCatch(solve(solution$hessian), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(hess_inv))
    warning("The Hessian at the solution is not positive definite (or could ",
            "not be inverted), so standard errors and confidence intervals ",
            "are NA.", call. = FALSE)
  # Standard errors on the internal scale.
  solution$par_stderr <- if (is.null(hess_inv)) {
    rep(NA_real_, length(solution$par))
  } else {
    sqrt(diag(hess_inv))
  }
  names(solution$par_stderr) <- par_names
  # 95% Wald confidence limits (estimate +/- 1.96 * SE) on the internal scale.
  # The steps below convert estimates, SEs and limits to the reporting scale.
  solution$par_upper <- solution$par + 1.96 * solution$par_stderr
  solution$par_lower <- solution$par - 1.96 * solution$par_stderr

  # log(alpha) on the original covariate scale. Centering moves a constant into
  # alpha, so log(alpha) = log(alpha') + sum((ctr / scl) * beta'), where alpha'
  # and beta' come from the standardized fit. This is a linear combination of the
  # fitted parameters, so its variance is exactly a' %*% vcov %*% a.
  a_vec <- c(0, 1, ctr / scl)
  log_alpha <- sum(a_vec * solution$par)
  se_log_alpha <- if (is.null(hess_inv)) {
    NA_real_
  } else {
    sqrt(drop(t(a_vec) %*% hess_inv %*% a_vec))
  }
  
  # Coefficients: convert from standardized units back to original covariate
  # units (centering does not change them; scaling by scl does)
  b_idx <- seq_along(solution$par)[-(1:2)]
  solution$par[b_idx]        <- solution$par[b_idx]        / scl
  solution$par_stderr[b_idx] <- solution$par_stderr[b_idx] / scl
  solution$par_upper[b_idx]  <- solution$par_upper[b_idx]  / scl
  solution$par_lower[b_idx]  <- solution$par_lower[b_idx]  / scl

  # r: exponentiate, then delta-method SE on the original scale,
  # SE(r) is approximately r * SE(log r)
  solution$par[1]        <- exp(solution$par[1])
  solution$par_upper[1]  <- exp(solution$par_upper[1])
  solution$par_lower[1]  <- exp(solution$par_lower[1])
  solution$par_stderr[1] <- solution$par[1] * solution$par_stderr[1]

  # alpha: original covariate scale, with the interval built on the log scale
  solution$par[2]        <- exp(log_alpha)
  solution$par_upper[2]  <- exp(log_alpha + 1.96 * se_log_alpha)
  solution$par_lower[2]  <- exp(log_alpha - 1.96 * se_log_alpha)
  solution$par_stderr[2] <- solution$par[2] * se_log_alpha

  # Keep the standardization constants (the hessian is on the standardized scale)
  solution$center <- ctr
  solution$scale  <- scl
  # Keep the terms (response removed) and factor levels so new data can be coded
  # the same way at prediction time (used by G2G_static_survival).
  tt <- terms(mf)
  solution$terms   <- delete.response(tt)
  solution$xlevels <- .getXlevels(tt, mf)
  # Tag the result with class "G2Gcov" (print method in methods.R).
  class(solution) <- c("G2Gcov", "list")
  return(solution)
}


#' Survival Curve for the G2G Model with Time-Invariant Covariates
#'
#' Computes the survival curve implied by a fitted static G2G model for one set
#' of covariate values:
#' \eqn{S(t) = (1 + t \cdot \exp(\mathbf{x}'\boldsymbol{\beta}) /
#' \alpha)^{-r}}, with \eqn{S(0) = 1}.
#'
#' @param solution A fitted model of class \code{"G2Gcov"} returned by
#'   \code{\link{G2G_static_MLE}}. Its \code{par} holds r, alpha and the
#'   coefficients on the original covariate scale (alpha already absorbs the
#'   centering used in fitting), so \code{X_vect} needs no centering or scaling.
#' @param X_vect Covariate values for one subject. Only one set of values is
#'   accepted per call; to get curves for several subjects or covariate
#'   profiles, call the function once per profile (for example in a loop or
#'   with \code{sapply}). Either a one-row data frame containing the original
#'   covariate variables (factors and dummy variables are coded as in the fit,
#'   using the stored terms and factor levels; a missing value is an error), or
#'   plain numbers: one value per model coefficient, as a numeric vector or
#'   one-row matrix named like the coefficients (any order) or unnamed in the
#'   same order as the coefficients. If the model has a single covariate
#'   coefficient, a single number. Anything else (more than one row, or the
#'   wrong number of values) is an error.
#' @param timeperiod Number of periods, a whole number >= 1. The curve is
#'   returned for periods 0, 1, ..., \code{timeperiod}.
#' @return Numeric vector of length \code{timeperiod + 1} with the survival
#'   probabilities S(0) = 1, S(1), ..., S(timeperiod).
#' @export
#' @importFrom stats model.frame model.matrix na.fail
#' @examples
#' \donttest{
#'   if (requireNamespace("survival", quietly = TRUE)) {
#'     library(survival)
#'     fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
#'     G2G_static_survival(fit, data.frame(age = 60, karno = 70),
#'                         timeperiod = 30)
#'   }
#' }
G2G_static_survival <- function(solution, X_vect, timeperiod) {
  # timeperiod sets how many periods are computed; reject values that
  # would silently return only S(0) = 1 or fail obscurely.
  if (!is.numeric(timeperiod) || length(timeperiod) != 1L || is.na(timeperiod) ||
      timeperiod < 1 || timeperiod != round(timeperiod))
    stop("'timeperiod' must be a single whole number >= 1.", call. = FALSE)

  # Fitted parameters (original covariate scale): r, alpha and the coefficients.
  r <- solution$par[1]
  alpha <- solution$par[2]
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

  # Periods 1, 2, ..., timeperiod.
  cen    <- seq_len(timeperiod)
  # Constant covariate effect exp(x' * beta) (the covariates do not change over time).
  C_c    <- c(exp(X %*% beta))
  # log S(t) = -r * log(1 + t * C / alpha) for each period t.
  LL_cen <- -r * log(1 + C_c * cen / alpha)

  # Prepend S(0) = 1 and exponentiate back to the survival scale.
  return(c(1, exp(LL_cen)))
}

#' Retention Rates for the G2G Model with Time-Invariant Covariates
#'
#' Computes the period-by-period retention rate implied by a fitted static G2G
#' model for one set of covariate values: the probability of surviving period t
#' given survival through period t - 1, \eqn{S(t) / S(t - 1)}, where \eqn{S} is
#' the curve from \code{\link{G2G_static_survival}}.
#'
#' @inheritParams G2G_static_survival
#' @return Numeric vector of length \code{timeperiod + 1}: a leading 1 (for
#'   time 0), then S(1)/S(0), S(2)/S(1), ..., S(T)/S(T-1) with
#'   T = \code{timeperiod}.
#' @export
#' @examples
#' \donttest{
#'   if (requireNamespace("survival", quietly = TRUE)) {
#'     library(survival)
#'     fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
#'     G2G_static_retention(fit, data.frame(age = 60, karno = 70),
#'                          timeperiod = 30)
#'   }
#' }
G2G_static_retention <- function(solution, X_vect, timeperiod) {
  # Survival curve c(S(0), ..., S(T)).
  st<-G2G_static_survival(solution, X_vect, timeperiod)
  # Ratio of consecutive survival values, with a leading 1 so the result lines
  # up with the survival vector.
  return(c(1, st[-1] / st[-length(st)]))
}
