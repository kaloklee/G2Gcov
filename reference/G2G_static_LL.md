# Log-Likelihood for G2G Model with Static Covariates

Internal function that computes the negative log-likelihood for the G2G
model with time-invariant covariates (one row per subject). Used by the
optimization routine in
[`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md).

## Usage

``` r
G2G_static_LL(par, df)
```

## Arguments

- par:

  Numeric vector of parameters on the optimizer's internal scale. The
  first two elements are log(r) and log(alpha) (both are exponentiated
  inside the function), and the remaining elements are the covariate
  coefficients.

- df:

  Data frame with columns, in this order: time (the observed period, a
  whole number \>= 1), status (1 = event, 0 = censored), then the
  covariates. One row per subject.

## Value

Scalar value of the negative log-likelihood
