# Conditional Survival Over Planned Future Periods Given History

For a customer who has survived the periods in a covariate history,
computes the probability of surviving each of the next `s` periods under
a planned covariate path: \\P(T \> t + k \mid T \> t)\\ for \\k = 1,
\dots, s\\.

## Usage

``` r
G2G_varying_cond_survival(solution, X_hist, X_future)
```

## Arguments

- solution:

  A fitted model of class `"G2Gcov"` returned by
  [`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md).
  Its `par` holds r, alpha and the coefficients on the original
  covariate scale (alpha already absorbs the centering used in fitting),
  so `X_mat` needs no centering or scaling.

- X_hist:

  Covariate values for each period the customer has survived, one row
  per period, in any of the forms accepted by `X_mat` in
  [`G2G_varying_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_survival.md).
  The history must contain only event-free periods. `NULL` (no history)
  gives the unconditional survival curve for the planned path.

- X_future:

  Planned covariate values for the periods to forecast, one row per
  future period (periods t + 1, ..., t + s), in the same form as
  `X_hist`. At least one row is required.

## Value

Numeric vector of length s + 1: a leading 1 (no additional periods),
then the probability of surviving 1, 2, ..., s additional periods given
survival through the history.

## Details

The customer's history updates the distribution of the baseline rate to
\\\mathrm{Gamma}(r, \alpha + C(t))\\, and averaging the future survival
over that posterior gives \$\$P(T \> t + k \mid T \> t) =
\left(\frac{\alpha + C(t)}{\alpha + C(t + k)}\right)^r,\$\$ where \\C(m)
= \sum\_{j=1}^{m} e^{z(j)\beta'}\\ accumulates the history and then the
planned covariates. This equals \\S(t + k) / S(t)\\ from
[`G2G_varying_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_survival.md)
for the combined path. The history length \\t\\ and the horizon \\s\\
are taken from the number of rows of `X_hist` and `X_future`.

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
  # Survived 8 weeks; plan 4 more weeks with a coupon and no other promotion
  hist   <- kb_data[kb_data$id == kb_data$id[1], c("coupon", "anyp")][1:8, ]
  future <- data.frame(coupon = rep(1, 4), anyp = rep(0, 4))
  G2G_varying_cond_survival(fit, hist, future)
#> Error in na.fail.default(structure(list(coupon = c(0, NA, NA, NA, NA, NA, NA, NA), anyp = c(55.61, NA, NA, NA, NA, NA, NA, NA)), class = "data.frame", row.names = c("1", "NA", "NA.1", "NA.2", "NA.3", "NA.4", "NA.5", "NA.6"), terms = ~coupon +     anyp)): missing values in object
# }
```
