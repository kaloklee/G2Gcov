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
  covariates. Transformations (such as `I(x - 60)`, `log(x)`, or
  `factor(x)`) and interactions are supported on the right-hand side,
  but the intercept must be kept (no `- 1` or `+ 0`): the intercept
  column is always dropped because \\\alpha\\ plays its role.

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
  data(kb_data_sample)
  fit <- G2G_varying_MLE(
    Surv(week, status) ~ coupon + anyp,
    data    = kb_data_sample,
    subject = "id"
  )
  print(fit)
#> G2G Model Fit
#> ---------------------------------------------------- 
#> Convergence:       Converged (code 0)
#> Log-likelihood:    -11.1869
#> AIC:               30.3739
#> Parameters:        4
#> 
#> Coefficients:
#>              Estimate Std.Error Lower.95   Upper.95
#> r (shape)      0.0611    0.1047   0.0021     1.7580
#> alpha (rate)   6.3229   31.1821   0.0004 99709.9104
#> coupon        24.9515   22.4348 -19.0199    68.9230
#> anyp          -0.0159    0.0554  -0.1246     0.0927
#> 
#> Note: Wide confidence intervals for r and alpha reflect a flat likelihood ridge (common with small samples, few observed events, or very low unobserved heterogeneity). Interpret r and alpha jointly (e.g., via plot_p0() or plot_lambda0()).

  # 90% confidence intervals instead of the default 95%
  fit90 <- G2G_varying_MLE(
    Surv(week, status) ~ coupon + anyp,
    data    = kb_data_sample,
    subject = "id",
    level   = 0.90
  )
# }
```
