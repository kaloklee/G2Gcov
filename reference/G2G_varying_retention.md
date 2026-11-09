# Retention Rates for the G2G Model with Time-Varying Covariates

Computes the period-by-period retention rate implied by a fitted G2G
model for one covariate path: the probability of surviving period t
given survival through period t - 1, `S(t) / S(t - 1)`, where `S` is the
curve from
[`G2G_varying_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_survival.md).

## Usage

``` r
G2G_varying_retention(solution, X_mat)
```

## Arguments

- solution:

  A fitted model of class `"G2Gcov"` returned by
  [`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md).
  Its `par` holds r, alpha and the coefficients on the original
  covariate scale (alpha already absorbs the centering used in fitting),
  so `X_mat` needs no centering or scaling.

- X_mat:

  Covariate values for each period, one row per period (t = 1, ..., T).
  Either a data frame containing the original covariate variables
  (factors and dummy variables are coded as in the fit, using the stored
  terms and factor levels; missing values are an error), or a matrix
  whose columns are named like the model coefficients (any order) or, if
  unnamed, already in the same order as the coefficients.

## Value

Numeric vector of length T + 1: a leading 1 (for time 0), then
S(1)/S(0), S(2)/S(1), ..., S(T)/S(T-1).

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
  path <- kb_data[kb_data$id == kb_data$id[1], c("coupon", "anyp")]
  G2G_varying_retention(fit, path)
#>           r (shape) 
#> 1.0000000 0.9208866 
# }
```
