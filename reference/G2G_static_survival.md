# Survival Curve for the G2G Model with Time-Invariant Covariates

Computes the survival curve implied by a fitted static G2G model for one
set of covariate values: \\S(t) = (1 + t \cdot
\exp(\mathbf{x}'\boldsymbol{\beta}) / \alpha)^{-r}\\, with \\S(0) = 1\\.

## Usage

``` r
G2G_static_survival(solution, X_vect, timeperiod)
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

Numeric vector of length `timeperiod + 1` with the survival
probabilities S(0) = 1, S(1), ..., S(timeperiod).

## Examples

``` r
# \donttest{
  if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
    G2G_static_survival(fit, data.frame(age = 60, karno = 70),
                        timeperiod = 30)
  }
#>  [1] 1.0000000 0.9928962 0.9858521 0.9788670 0.9719403 0.9650716 0.9582602
#>  [8] 0.9515057 0.9448073 0.9381647 0.9315773 0.9250445 0.9185658 0.9121407
#> [15] 0.9057687 0.8994494 0.8931821 0.8869664 0.8808018 0.8746878 0.8686240
#> [22] 0.8626098 0.8566449 0.8507286 0.8448607 0.8390405 0.8332678 0.8275419
#> [29] 0.8218626 0.8162293 0.8106416
# }
```
