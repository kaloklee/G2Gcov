# Fit a G2G Model with Time-Invariant Covariates

Fits a Grassia(II)-Geometric (G2G) survival model with time-invariant
covariates using maximum likelihood estimation via BFGS. Each subject is
represented by a single row. The survivor function at discrete period
\\t\\ is \\S(t) = (1 + t \cdot \exp(\mathbf{x}'\boldsymbol{\beta}) /
\alpha)^{-r}\\, where \\r\\ and \\\alpha\\ are the shape and rate of the
Gamma mixing distribution and \\\boldsymbol{\beta}\\ are covariate
log-hazard coefficients.

Covariates are centered and scaled internally to keep the likelihood
numerically stable, and the results are converted back to the original
covariate scale.

## Usage

``` r
G2G_static_MLE(formula, data, level = 0.95)
```

## Arguments

- formula:

  A formula of the form `Surv(time, status) ~ x1 + x2 + ...` where
  `time` is the observed period, `status` is the event indicator (1 =
  event, 0 = censored), and `x1`, `x2`, etc. are time-invariant
  covariates. Transformations, interactions and factors are allowed on
  the right-hand side, but the intercept must be kept (no `- 1` or
  `+ 0`): the intercept column is always dropped, because \\\alpha\\
  plays its role.

- data:

  A data frame containing all variables referenced in `formula`. One row
  per subject. The observed times must be whole numbers \>= 1 (discrete
  periods).

- level:

  Confidence level for the Wald confidence intervals, a single number
  strictly between 0 and 1. Default is 0.95.

## Value

A list of class `"G2Gcov"` containing the optimization results from
[`optim`](https://rdrr.io/r/stats/optim.html) with the following
additional elements. All estimates, standard errors and confidence
bounds are reported on the original covariate scale (only `hessian` is
on the internal parameterization):

- par:

  Named numeric vector of parameter estimates. The first two elements
  are `r (shape)` and `alpha (rate)`, returned on the natural
  (exponentiated) scale. Remaining elements are covariate log-hazard
  coefficients in original units.

- par_stderr:

  Standard errors on the original scale (named like `par`). For `r` and
  `alpha` these come from the delta method (e.g. SE(r) is approximately
  r \* SE(log r)). If the Hessian at the optimum is not positive
  definite (or cannot be inverted) a warning is issued and these are
  `NA`.

- par_upper:

  Upper bounds of the Wald confidence intervals at the requested
  `level`. For `r` and `alpha` the interval is formed on the log scale
  and exponentiated.

- par_lower:

  Lower bounds of the Wald confidence intervals at the requested
  `level`.

- level:

  The confidence level used for `par_lower` and `par_upper`.

- value:

  Negative log-likelihood at the optimum.

- convergence:

  Convergence code from `optim` (0 = success; a warning is issued for
  any other code).

- hessian:

  Observed Hessian matrix at the optimum, on the internal
  parameterization (log r, log alpha, standardized coefficients), with
  row and column names.

- center:

  Covariate means used internally for centering.

- scale:

  Covariate standard deviations used internally for scaling.

- terms:

  The `terms` object for the covariates (response removed), used to code
  new data in
  [`G2G_static_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_survival.md).

- xlevels:

  Factor levels of the covariates at fitting time, also used to code new
  data.

## See also

[`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md)
for the time-varying extension;
[`G2G_static_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_survival.md)
and
[`G2G_static_retention`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_retention.md)
for the fitted curves;
[`plot_p0`](https://kaloklee.github.io/G2Gcov/reference/plot_p0.md) to
visualize the fitted baseline propensity distribution.

## Examples

``` r
# \donttest{
  if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
    print(fit)

    # 90% confidence intervals instead of the default 95%
    fit90 <- G2G_static_MLE(
      Surv(time, status) ~ age + karno,
      data  = veteran,
      level = 0.90
    )
  }
#> G2G Model Fit
#> ---------------------------------------------------- 
#> Convergence:       Converged (code 0)
#> Log-likelihood:    -723.3544
#> AIC:               1454.7087
#> Parameters:        4
#> 
#> Coefficients:
#>              Estimate Std.Error Lower.95 Upper.95
#> r (shape)      5.4654    3.2150   1.7254  17.3117
#> alpha (rate)  48.0151   52.4109   5.6526 407.8554
#> age           -0.0027    0.0102  -0.0228   0.0174
#> karno         -0.0373    0.0054  -0.0478  -0.0267
# }
```
