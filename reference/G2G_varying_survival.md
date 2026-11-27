# Survival Curve for the G2G Model with Time-Varying Covariates

Computes the survival curve implied by a fitted G2G model for one
covariate path: `S(t) = (1 + C(t) / alpha)^(-r)`, where `C(t)` is the
cumulative sum, over periods 1 to t, of the exponentiated linear
predictor (covariates times coefficients), and `S(0) = 1`.

## Usage

``` r
G2G_varying_survival(solution, X_mat)
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

Numeric vector of length T + 1 with the survival probabilities S(0) = 1,
S(1), ..., S(T).

## Examples

``` r
# \donttest{
  data(kb_data_sample)
  fit <- G2G_varying_MLE(
    Surv(week, status) ~ coupon + anyp,
    data    = kb_data_sample,
    subject = "id"
  )
  # Survival curve for the covariate path of Subject 165
  path <- kb_data_sample[kb_data_sample$id == 165, c("coupon", "anyp")]
  G2G_varying_survival(fit, path)
#>  [1] 1.0000000 0.9961475 0.9925494 0.9873133 0.9517759 0.9435223 0.9403493
#>  [8] 0.9377270 0.9339596 0.9329796 0.9316340 0.9229708 0.9181504 0.9138431
#> [15] 0.9113546 0.9102154 0.9092013 0.9082543 0.9073679 0.9061442 0.9046525
#> [22] 0.9029872 0.9013737 0.9004307 0.8995220
# }
```
