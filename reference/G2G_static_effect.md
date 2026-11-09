# Covariate Effect for One Subject with Time-Invariant Covariates

Internal helper shared by
[`G2G_static_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_survival.md),
[`G2G_static_posterior_p`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_posterior_p.md)
and
[`G2G_static_cond_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_cond_survival.md).
Codes one set of covariate values using the terms and factor levels
stored in the fit and returns the constant covariate effect
\\\exp(x'\beta)\\.

## Usage

``` r
G2G_static_effect(solution, X_vect)
```

## Arguments

- solution:

  A fitted model of class `"G2Gcov"` returned by
  [`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md).

- X_vect:

  Covariate values for one subject, in any of the forms accepted by
  `X_vect` in
  [`G2G_static_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_survival.md).
  Only one set of values is accepted; anything else is an error.

## Value

A single number, \\\exp(x'\beta)\\.
