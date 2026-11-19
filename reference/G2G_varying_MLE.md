# Fit G2G Model with Time-Varying Covariates

Fits a Grassia(II)-Geometric (G2G) survival model with time-varying
covariates using maximum likelihood estimation. Covariates are centered
and scaled internally to keep the likelihood numerically stable, and the
results are converted back to the original covariate scale.

## Usage

``` r
G2G_varying_MLE(fo, data, subject, level = 0.95)
```

## Arguments

- fo:

  Formula object specifying the model in the format
  `Surv(time, status) ~ x1 + x2 + ...` where time is the observation
  time, status is the event indicator, and x1, x2, etc. are time-varying
  covariates. Note: The formula should use variable names, not column
  indices. Only plain variable names are allowed on the right-hand side:
  transformations such as `log(x)`, `I(x^2)` or `factor(x)` and
  interactions such as `x1 * x2` are rejected with an error (they would
  otherwise be silently ignored). Create transformed variables as new
  columns in `data` first.

- data:

  Data frame containing all variables specified in the formula. Data
  should be in long format with one row per time point per subject
  (person-period format). Rows are sorted by subject and time
  internally. Within each subject, the time variable must run 1, 2, 3,
  ... with no gaps or duplicates. A subject counts as having had an
  event if the status variable is 1 on any of their rows, and the event
  is placed in the subject's last period.

- subject:

  Character string specifying the name of the subject ID column in the
  data frame. This column identifies which observations belong to the
  same subject over time.

- level:

  Confidence level for the Wald confidence intervals, a single number
  strictly between 0 and 1. Default is 0.95.

## Value

A list of class `"G2Gcov"` containing the optimization results from
[`optim`](https://rdrr.io/r/stats/optim.html) with additional elements.
All estimates, standard errors and confidence bounds are reported on the
original covariate scale (only `hessian` is on the internal
parameterization):

- par:

  Parameter estimates. First two are r and alpha (shape and rate
  parameters), remaining are covariate coefficients in original units.

- par_stderr:

  Standard errors of parameter estimates (named like `par`). For r and
  alpha these come from the delta method, e.g. SE(r) is approximately r
  \* SE(log r). If the Hessian is not positive definite (or cannot be
  inverted) a warning is issued and these are `NA`.

- par_upper:

  Upper bounds of the Wald confidence intervals at the requested
  `level`. For r and alpha the interval is formed on the log scale and
  exponentiated.

- par_lower:

  Lower bounds of the Wald confidence intervals at the requested `level`

- level:

  The confidence level used for `par_lower` and `par_upper`

- value:

  Negative log-likelihood at the optimum

- convergence:

  Convergence code (0 indicates successful convergence; a warning is
  issued for any other code)

- hessian:

  Hessian matrix at the optimum, on the internal parameterization (log
  r, log alpha, standardized coefficients), with row and column names

- center:

  Covariate means used internally for centering

- scale:

  Covariate standard deviations used internally for scaling

- terms:

  The `terms` object for the covariates, used to code new data in
  [`G2G_varying_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_survival.md)

- xlevels:

  Factor levels of the covariates at fitting time, also used to code new
  data

## Examples

``` r
# \donttest{
  data(kb_data)
  kb_data$status <- 1L - kb_data$censor
  fit <- G2G_varying_MLE(
    "Surv(week, status) ~ coupon + anyp",
    data    = kb_data,
    subject = "id"
  )
#> Warning: The Hessian at the solution is not positive definite (or could not be inverted), so standard errors and confidence intervals are NA.
  print(fit)
#> G2G Model Fit
#> ---------------------------------------------------- 
#> Convergence:       Converged (code 0)
#> Log-likelihood:    -5612.1808
#> AIC:               11232.3616
#> Parameters:        4
#> 
#> Coefficients:
#>                Estimate Std.Error Lower.95 Upper.95
#> r (shape)     2996.6875        NA       NA       NA
#> alpha (rate) 40722.0679        NA       NA       NA
#> coupon         -78.2544        NA       NA       NA
#> anyp             0.0020        NA       NA       NA
#> 
#> Note: some standard errors are NA (the Hessian at the optimum is singular or not positive definite).

  # 90% confidence intervals instead of the default 95%
  fit90 <- G2G_varying_MLE(
    "Surv(week, status) ~ coupon + anyp",
    data    = kb_data,
    subject = "id",
    level   = 0.90
  )
#> Warning: The Hessian at the solution is not positive definite (or could not be inverted), so standard errors and confidence intervals are NA.
# }
```
