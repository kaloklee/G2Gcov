#' Log-Likelihood for G2G Model with Time-Varying Covariates
#'
#' Internal function that computes the negative log-likelihood for the G2G model
#' with time-varying covariates. Used by the optimization routine.
#'
#' @param par Numeric vector of parameters on the optimizer's internal scale.
#'   The first two elements are log(r) and log(alpha) (the logs of the shape and
#'   rate parameters; both are exponentiated inside the function), and the
#'   remaining elements are the covariate coefficients.
#' @param data_df Data frame with columns, in this order:
#'   \describe{
#'     \item{id}{Subject identifier}
#'     \item{time}{Time point for each observation. Within each subject, rows
#'       must be ordered by time, with time running 1, 2, 3, ... and no gaps}
#'     \item{status}{Event indicator (0 = no occurrence, 1 = event occurred).
#'       A subject counts as having had an event if status is 1 on any of
#'       their rows, and the event is placed in the subject's last period}
#'     \item{...}{Time-varying covariates (remaining columns)}
#'   }
#' @return Scalar value of the negative log-likelihood
#' @keywords internal
G2G_varying_LL <- function(par,data_df) {
  # Overview: the survival function implied by the G2G model is
  #   S(t) = (1 + C(t) / alpha)^(-r),   with S(0) = 1,
  # where C(t) = sum over periods 1..t of exp(x_s' * coeff) is the cumulative
  # covariate effect for a subject. Each subject contributes to the likelihood:
  #   - event in period t:   log[ S(t-1) - S(t) ]
  #   - censored at time t:  log[ S(t) ]
  # The function returns the NEGATIVE of the summed log-likelihood so that
  # optim() can minimize it.
  #
  # Expected layout of data_df: columns 1-3 are id, time and status; columns 4
  # onward are the covariates. One row per subject per period, with time
  # running 1, 2, 3, ... for each subject.
  #
  # par = c(log r, log alpha, covariate coefficients). The commented-out lines
  # below record an earlier mean/polarization parameterization.
  #par[1] = mean of BG
  #par[2] = polarization of BG
  #r=par[1]*(1/par[2]-1); 
  #alpha=(1-par[1])*(1/par[2]-1);
  
  # Back-transform r and alpha from the log scale used by the optimizer
  # (guarantees both are positive); the rest of par are the coefficients.
  r = exp(par[1]);
  alpha = exp(par[2]);
  coeff=par[-(1:2)];
  
  #X = scale(as.matrix(data_df[,-(1:3)]));
  # Covariate matrix: everything after the id, time and status columns.
  X = as.matrix(data_df[,-(1:3)]);
  
  #print(length(coeff)==dim(X)[2])
  
  # Per-period covariate effect: exp(x_t' * coeff).
  data_df$Ct = exp(X %*% coeff);
  
  # Running total of the effect within each subject, i.e. C(t) for every
  # person-period row.
  data_df$cumsumCt <- ave(data_df$Ct, data_df$id, FUN = cumsum);
  
  #gather unique id information
  
#  print(names(data_df))
  
  # One row per subject: last observed time (max of time) and event
  # indicator (max of status, which is 1 if the subject's event occurred).
  id_df <- aggregate(cbind(time, status) ~ id, data = data_df, FUN = max)
  
  #uncensored piece of likelihood
  
  # Running total of the log-likelihood from subjects with an observed event.
  LL_uncen = 0;
  
  # Distinct times at which at least one subject had an event. Subjects are
  # processed in groups by event time so the work can be vectorized.
  uncen_time = unique(id_df[which(id_df$status==1),"time"]);
  
  for (t in uncen_time) {
    
    # Subjects whose event occurred at time t.
    id_need = id_df[which(id_df$status == 1 & id_df$time == t),"id"];
    
    # Event in the first period: S(0) = 1, so P(event) = 1 - S(1). Otherwise:
    # P(event in period t) = S(t-1) - S(t).
    if (t == 1) {
      
      # C(1) for these subjects.
      cumsumCt_b = data_df[which(data_df$id %in% id_need & data_df$time ==1),
                           "cumsumCt"]; 
      
      # log(1 - S(1)), written as logdiffexp(0, log S(1)). logdiffexp(a, b) is
      # log(exp(a) - exp(b)) computed stably on the log scale (defined in utils.R).
      LL_uncen = LL_uncen+sum(logdiffexp(0,-r*log(1+cumsumCt_b/alpha)));
      
    }
    
    else { 
      
      # C(t-1) (cumsumCt_a) and C(t) (cumsumCt_b) for these subjects.
      cumsumCt_a = data_df[which(data_df$id %in% id_need & data_df$time == (t-1)),
                           "cumsumCt"]; 
      cumsumCt_b = data_df[which(data_df$id %in% id_need & data_df$time == (t)),
                           "cumsumCt"]; 
      
      
      # log(S(t-1) - S(t)) for each subject in the group, summed into the total.
      LL_uncen= LL_uncen + sum(logdiffexp(-r*log(1+cumsumCt_a/alpha), 
                                      -r*log(1+cumsumCt_b/alpha)));  
      
    }
    
  }
  
  #censored piece of likelihood  
  
  # Distinct times at which at least one subject was censored (no event seen).
  cen_time = unique(id_df[which(id_df$status==0),"time"]);
  
  # Running total of the log-likelihood from censored subjects.
  LL_cen = 0;
  
  for (t in cen_time) {
    
    # Subjects censored at time t.
    id_need = id_df[which(id_df$status == 0 & id_df$time == t),"id"];
    
    # C(t) for these subjects.
    cumsumCt_b = data_df[which(data_df$id %in% id_need & data_df$time == (t)),
                         "cumsumCt"]; 
    
    # A censored subject survived through period t, so the contribution is
    # log S(t) = -r * log(1 + C(t) / alpha).
    LL_cen= LL_cen + sum((-r)*log(1+cumsumCt_b/alpha));
    
  }  
  
  #  print ( -(LL_uncen) )
  #  print( -(LL_cen) )
  #  print(par)
  
  # Negative log-likelihood (optim() minimizes, so the sign is flipped).
  return ( -(LL_uncen)-(LL_cen) );
  
}

#' Fit G2G Model with Time-Varying Covariates
#'
#' Fits a Grassia(II)-Geometric (G2G) survival model with time-varying covariates
#' using maximum likelihood estimation. Covariates are centered and scaled
#' internally to keep the likelihood numerically stable, and the results are
#' converted back to the original covariate scale.
#'
#' @param fo Formula object specifying the model in the format 
#'   \code{Surv(time, status) ~ x1 + x2 + ...} where time is the observation
#'   time, status is the event indicator, and x1, x2, etc. are time-varying
#'   covariates. Note: The formula should use variable names, not column indices.
#'   Only plain variable names are allowed on the right-hand side: transformations
#'   such as \code{log(x)}, \code{I(x^2)} or \code{factor(x)} and interactions
#'   such as \code{x1 * x2} are rejected with an error (they would otherwise be
#'   silently ignored). Create transformed variables as new columns in
#'   \code{data} first.
#' @param data Data frame containing all variables specified in the formula.
#'   Data should be in long format with one row per time point per subject
#'   (person-period format). Rows are sorted by subject and time internally.
#'   Within each subject, the time variable must run 1, 2, 3, ... with no gaps
#'   or duplicates. A subject counts as having had an event if the status
#'   variable is 1 on any of their rows, and the event is placed in the
#'   subject's last period.
#' @param subject Character string specifying the name of the subject ID column
#'   in the data frame. This column identifies which observations belong to
#'   the same subject over time.
#'
#' @return A list of class \code{"G2Gcov"} containing the optimization results
#'   from \code{\link[stats]{optim}} with additional elements. All estimates,
#'   standard errors and confidence bounds are reported on the original covariate
#'   scale (only \code{hessian} is on the internal parameterization):
#'   \item{par}{Parameter estimates. First two are r and alpha (shape and rate
#'     parameters), remaining are covariate coefficients in original units.}
#'   \item{par_stderr}{Standard errors of parameter estimates (named like
#'     \code{par}). For r and alpha these come from the delta method, e.g. SE(r)
#'     is approximately r * SE(log r). If the Hessian is not positive definite
#'     (or cannot be inverted) a warning is issued and these are \code{NA}.}
#'   \item{par_upper}{Upper bounds of 95\% confidence intervals. For r and alpha
#'     the interval is formed on the log scale and exponentiated.}
#'   \item{par_lower}{Lower bounds of 95\% confidence intervals}
#'   \item{value}{Negative log-likelihood at the optimum}
#'   \item{convergence}{Convergence code (0 indicates successful convergence; a
#'     warning is issued for any other code)}
#'   \item{hessian}{Hessian matrix at the optimum, on the internal
#'     parameterization (log r, log alpha, standardized coefficients), with
#'     row and column names}
#'   \item{center}{Covariate means used internally for centering}
#'   \item{scale}{Covariate standard deviations used internally for scaling}
#'   \item{terms}{The \code{terms} object for the covariates, used to code new
#'     data in \code{\link{G2G_varying_survival}}}
#'   \item{xlevels}{Factor levels of the covariates at fitting time, also used
#'     to code new data}
#'
#' @export
#' @importFrom stats aggregate optim as.formula ave model.matrix model.frame terms
#' @importFrom stats .getXlevels
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE(
#'     "Surv(week, status) ~ coupon + anyp",
#'     data    = kb_data,
#'     subject = "id"
#'   )
#'   print(fit)
#' }
#'
G2G_varying_MLE <- function(fo, data, subject) {

  # Workflow: (1) validate the inputs, (2) build and standardize the covariate
  # matrix, (3) maximize the likelihood via G2G_varying_optim(), (4) get standard
  # errors and 95% Wald intervals from the inverse Hessian, and (5) convert the
  # results back to the original covariate scale.
  #
  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(subject) || length(subject) != 1L)
    stop("'subject' must be a single character string.", call. = FALSE)
  if (!subject %in% names(data))
    stop("Subject column '", subject, "' not found in 'data'.", call. = FALSE)

  # Variable names from the formula, in order: time, status, then covariates
  # (for Surv(time, status) ~ x1 + x2 this is time, status, x1, x2).
  vars        <- all.vars(as.formula(fo))
  time_name   <- vars[1]
  status_name <- vars[2]
  cov_names   <- vars[-(1:2)]

  if (length(cov_names) == 0L)
    stop("Formula must include at least one covariate on the right-hand side.",
         call. = FALSE)

  missing_cols <- setdiff(c(time_name, status_name, cov_names), names(data))
  if (length(missing_cols) > 0L)
    stop("Column(s) not found in 'data': ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)

  # The formula is only used to read off variable names, so every right-hand-side
  # term must be a plain variable. Transformations (log(x), I(x^2), factor(x))
  # and interactions (x1 * x2) would otherwise be silently ignored.
  term_labels <- gsub("`", "", attr(terms(as.formula(fo)), "term.labels"),
                      fixed = TRUE)
  if (!setequal(term_labels, cov_names))
    stop("The right-hand side of the formula must contain only plain variable ",
         "names (no transformations or interactions). Offending term(s): ",
         paste(setdiff(union(term_labels, cov_names),
                       intersect(term_labels, cov_names)), collapse = ", "),
         ". Create transformed variables as new columns in 'data' first.",
         call. = FALSE)

  na_cols <- c(subject, time_name, status_name, cov_names)
  na_cols <- na_cols[vapply(na_cols,
                            function(v) anyNA(data[[v]]), logical(1L))]
  if (length(na_cols) > 0L)
    stop("NA values found in column(s): ",
         paste(na_cols, collapse = ", "),
         ". Remove or impute missing values before fitting.", call. = FALSE)

  if (!all(unique(data[[status_name]]) %in% c(0, 1)))
    stop("Status column '", status_name,
         "' must contain only 0 and 1.", call. = FALSE)
  if (sum(data[[status_name]]) == 0L)
    stop("No events found in '", status_name,
         "' (all values are 0). At least one event is required.", call. = FALSE)
  if (any(data[[time_name]] < 1))
    stop("Time column '", time_name,
         "' must contain values >= 1.", call. = FALSE)

  # Person-period layout. The likelihood needs each subject's rows in time order
  # with time running 1, 2, ..., n, because the cumulative sums follow row order.
  # Sorting here makes the input row order irrelevant; gaps and duplicated times
  # cannot be repaired automatically, so they are errors.
  data     <- data[order(data[[subject]], data[[time_name]]), , drop = FALSE]
  period   <- ave(seq_len(nrow(data)), data[[subject]], FUN = seq_along)
  bad_time <- data[[time_name]] != period
  if (any(bad_time)) {
    bad_ids <- unique(data[[subject]][bad_time])
    stop("Within each subject, time column '", time_name,
         "' must run 1, 2, 3, ... with no gaps or duplicates. Check subject(s): ",
         paste(bad_ids[seq_len(min(5L, length(bad_ids)))], collapse = ", "),
         ".", call. = FALSE)
  }
  # ---------------------------------------------------------------------------
  
  # Design matrix for the covariates only. model.matrix() expands factors into
  # dummy columns and adds an intercept column, which is dropped when
  # standardizing below (an intercept would be redundant with alpha).
  #df_temp<-data[,c(all.vars(as.formula(fo))[-c(1:2)])];
  df_temp <- data[, c(all.vars(as.formula(fo))[-c(1:2)]), drop = FALSE];
  X<-model.matrix( ~ ., data = df_temp);
  # Keep the terms and factor levels so new data can be coded the same way
  # at prediction time (used by G2G_varying_survival).
  mf_cov <- model.frame(~ ., data = df_temp)
  tt     <- attr(mf_cov, "terms")

  # Standardize the covariates (center and scale) so exp(X %*% coeff) cannot
  # overflow and the optimizer is well conditioned. This is an exact
  # reparameterization of the model: centering is absorbed by alpha and scaling
  # rescales the coefficients. Both are undone below. drop = FALSE keeps a matrix
  # when there is only one covariate column.
  Xs  <- scale(X[, -1, drop = FALSE])
  ctr <- attr(Xs, "scaled:center")
  scl <- attr(Xs, "scaled:scale")
  if (any(scl == 0))
    stop("Covariate(s) with no variation: ",
         paste(names(scl)[scl == 0], collapse = ", "),
         ". Remove them before fitting.", call. = FALSE)
  
  # Long-format data frame in the layout G2G_varying_LL() expects:
  # id, time, status, then the standardized covariates.
  model_data<-data.frame(id = unlist(data[,c(subject)]), 
                         time = unlist(data[,c(time_name)]),
                         status = unlist(data[,c(status_name)]),
                         X=Xs);
  # Maximize the likelihood. At this point the parameters are on the internal
  # scale: log r, log alpha and coefficients for the standardized covariates.
  solution <- G2G_varying_optim(model_data)

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
  # Covariance matrix of the internal-scale estimates: the inverse of the
  # Hessian of the negative log-likelihood. The inverse is only used when the
  # Hessian is finite and positive definite (all eigenvalues > 0), i.e. the
  # solution is a proper minimum. Otherwise the standard errors are returned as
  # NA with a warning.
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
  solution$terms   <- tt
  solution$xlevels <- .getXlevels(tt, mf_cov)
  # Tag the result with class "G2Gcov".
  class(solution) <- c("G2Gcov", "list")
  return(solution)

}

#' Internal Optimization for Time-Varying Model
#'
#' Internal wrapper function that performs the actual optimization for the
#' G2G model with time-varying covariates using the BFGS method. The shape and
#' rate parameters are optimized on the log scale, so an unconstrained method
#' is appropriate.
#'
#' @param model_data Data frame prepared by \code{G2G_varying_MLE} containing
#'   id, time, status, and the (standardized) covariate columns
#' @return The unmodified list returned by \code{\link[stats]{optim}}, with
#'   elements including:
#'   \item{par}{Parameter estimates on the internal scale: log(r), log(alpha)
#'     and the coefficients of the standardized covariates}
#'   \item{value}{Negative log-likelihood at the optimum}
#'   \item{convergence}{Convergence code (0 indicates successful convergence)}
#'   \item{hessian}{Hessian of the negative log-likelihood at the optimum}
#'   Standard errors and confidence intervals are not computed here;
#'   \code{G2G_varying_MLE} adds them.
#' @keywords internal
G2G_varying_optim <- function(model_data) {
  
  # Total number of parameters: log r, log alpha and one coefficient per
  # covariate (model_data has id, time, status plus the covariate columns).
  nvar=dim(model_data)[2]-1;
  
  # Starting values: r = alpha = 0.5 (passed to optim as log(0.5), because the
  # parameters are on the log scale) and all coefficients 0 (no covariate
  # effect). BFGS is an unconstrained method, which works here because
  # r and alpha are optimized on the log scale. hessian = TRUE returns the Hessian
  # of the negative log-likelihood, which the caller uses for standard errors.
  # The warning handler hides "NaNs produced" warnings, which typically come from
  # extreme trial parameter values during the search; other warnings pass through.
  solution <- withCallingHandlers(
    optim(par     = c(log(0.5), log(0.5), rep(0, nvar - 2)),
          fn      = G2G_varying_LL,
          data_df = model_data,
          method  = "BFGS",
          control = list(maxit = 1000),
          hessian = TRUE),
    warning = function(w) {
      if (grepl("NaNs produced", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )

  # Raw optim() output; the caller adds names, standard errors and intervals.
  return(solution)
  
}

#' Survival Curve for the G2G Model with Time-Varying Covariates
#'
#' Computes the survival curve implied by a fitted G2G model for one covariate
#' path: \code{S(t) = (1 + C(t) / alpha)^(-r)}, where \code{C(t)} is the
#' cumulative sum, over periods 1 to t, of the exponentiated linear predictor
#' (covariates times coefficients), and \code{S(0) = 1}.
#'
#' @param solution A fitted model of class \code{"G2Gcov"} returned by
#'   \code{\link{G2G_varying_MLE}}. Its \code{par} holds r, alpha and the
#'   coefficients on the original covariate scale (alpha already absorbs the
#'   centering used in fitting), so \code{X_mat} needs no centering or scaling.
#' @param X_mat Covariate values for each period, one row per period
#'   (t = 1, ..., T). Either a data frame containing the original covariate
#'   variables (factors and dummy variables are coded as in the fit, using the
#'   stored terms and factor levels; missing values are an error), or a matrix
#'   whose columns are named like the model coefficients (any order) or, if
#'   unnamed, already in the same order as the coefficients.
#' @return Numeric vector of length T + 1 with the survival probabilities
#'   S(0) = 1, S(1), ..., S(T).
#' @export
#' @importFrom stats model.frame model.matrix na.fail
#' @examples
#' \donttest{
#'   data(kb_data)
#'   kb_data$status <- 1L - kb_data$censor
#'   fit <- G2G_varying_MLE(
#'     "Surv(week, status) ~ coupon + anyp",
#'     data    = kb_data,
#'     subject = "id"
#'   )
#'   # Survival curve for the covariate path of the first subject
#'   path <- kb_data[kb_data$id == kb_data$id[1], c("coupon", "anyp")]
#'   G2G_varying_survival(fit, path)
#' }
#'
G2G_varying_survival <- function(solution, X_mat) {
  # Fitted parameters (original covariate scale).
  r     <- solution$par[1]
  alpha <- solution$par[2]
  beta  <- solution$par[-(1:2)]
  
  # Build the covariate matrix X (one row per period).
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
        stop("Column names of 'X_mat' must be: ",
             paste(names(beta), collapse = ", "), call. = FALSE)
      X <- X[, names(beta), drop = FALSE]      # reorder to match beta
    } else if (ncol(X) != length(beta)) {
      # Unnamed columns: only the number of columns can be checked, so they must
      # already be in the coefficient order.
      stop("'X_mat' must have ", length(beta), " columns.", call. = FALSE)
    }
  }
  
  # Per-period covariate effect exp(x_t' * beta); c() turns the one-column
  # matrix into a vector.
  Ct       <- c(exp(X %*% beta))
  # Cumulative effect C(t).
  cumsumCt <- cumsum(Ct)
  # S(0) = 1 followed by S(1..T) = exp(-r * log(1 + C(t) / alpha)).
  c(1, exp(-r * log(1 + cumsumCt / alpha)))
}

#' Retention Rates for the G2G Model with Time-Varying Covariates
#'
#' Computes the period-by-period retention rate implied by a fitted G2G model
#' for one covariate path: the probability of surviving period t given survival
#' through period t - 1, \code{S(t) / S(t - 1)}, where \code{S} is the curve
#' from \code{\link{G2G_varying_survival}}.
#'
#' @inheritParams G2G_varying_survival
#' @return Numeric vector of length T + 1: a leading 1 (for time 0), then
#'   S(1)/S(0), S(2)/S(1), ..., S(T)/S(T-1).
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
#'   path <- kb_data[kb_data$id == kb_data$id[1], c("coupon", "anyp")]
#'   G2G_varying_retention(fit, path)
#' }
#'
G2G_varying_retention <- function(solution, X_mat) {
  # Survival curve c(S(0), ..., S(T)).
  st<-G2G_varying_survival(solution, X_mat)
  # Ratio of consecutive survival values, with a leading 1 so the result lines
  # up with the survival vector.
  return(c(1, st[-1] / st[-length(st)]))
}

             
