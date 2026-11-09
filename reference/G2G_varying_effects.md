# Per-Period Covariate Effects for a Covariate Path

Internal helper shared by
[`G2G_varying_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_survival.md),
[`G2G_varying_posterior_p0`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_posterior_p0.md)
and
[`G2G_varying_cond_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_cond_survival.md).
Codes a covariate path using the terms and factor levels stored in the
fit and returns the per-period covariate effect \\\exp(x_j'\beta)\\, one
value per row of the path.

## Usage

``` r
G2G_varying_effects(solution, X_mat, arg = "X_mat")
```

## Arguments

- solution:

  A fitted model of class `"G2Gcov"` returned by
  [`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md).

- X_mat:

  Covariate path with one row per period: a data frame holding the
  original covariate variables (factors and dummies are coded as in the
  fit; missing values are an error), or a matrix whose columns are named
  like the model coefficients (any order) or, if unnamed, already in
  coefficient order. `NULL` or a path with zero rows means "no periods"
  and returns `numeric(0)`.

- arg:

  Name of the argument being checked, used in error messages.

## Value

Numeric vector with one effect per period (`numeric(0)` for an empty
path).
