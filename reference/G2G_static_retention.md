# Retention Rates for the G2G Model with Time-Invariant Covariates

Computes the period-by-period retention rate implied by a fitted static
G2G model for one set of covariate values: the probability of surviving
period t given survival through period t - 1, \\S(t) / S(t - 1)\\, where
\\S\\ is the curve from
[`G2G_static_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_survival.md).

## Usage

``` r
G2G_static_retention(solution, X_vect, timeperiod)
```

## Arguments

- solution:

  A fitted model of class `"G2Gcov"` returned by
  [`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md).
  Its `par` holds r, alpha and the coefficients on the original
  covariate scale (alpha already absorbs the centering used in fitting),
  so `X_vect` needs no centering or scaling.

- X_vect:

  Covariate values for one subject. Only one set of values is accepted
  per call; to get curves for several subjects or covariate profiles,
  call the function once per profile (for example in a loop or with
  `sapply`). Either a one-row data frame containing the original
  covariate variables (factors and dummy variables are coded as in the
  fit, using the stored terms and factor levels; a missing value is an
  error), or plain numbers: one value per model coefficient, as a
  numeric vector or one-row matrix named like the coefficients (any
  order) or unnamed in the same order as the coefficients. If the model
  has a single covariate coefficient, a single number. Anything else
  (more than one row, or the wrong number of values) is an error.

- timeperiod:

  Number of periods, a whole number \>= 1. The curve is returned for
  periods 0, 1, ..., `timeperiod`.

## Value

Numeric vector of length `timeperiod + 1`: a leading 1 (for time 0),
then S(1)/S(0), S(2)/S(1), ..., S(T)/S(T-1) with T = `timeperiod`.

## Examples

``` r
# \donttest{
  if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
    G2G_static_retention(fit, data.frame(age = 60, karno = 70),
                         timeperiod = 30)
  }
#>  [1] 1.0000000 0.9928962 0.9929055 0.9929147 0.9929238 0.9929330 0.9929421
#>  [8] 0.9929512 0.9929603 0.9929693 0.9929784 0.9929874 0.9929964 0.9930053
#> [15] 0.9930143 0.9930232 0.9930321 0.9930410 0.9930498 0.9930586 0.9930674
#> [22] 0.9930762 0.9930850 0.9930937 0.9931024 0.9931111 0.9931198 0.9931284
#> [29] 0.9931371 0.9931457 0.9931543
# }
```
