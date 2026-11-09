# Conditional Survival Over Future Periods Given History (Time-Invariant Covariates)

For a subject who has survived `hist` periods, computes the probability
of surviving each of the next `future` periods: \\P(T \> \mathrm{hist} +
k \mid T \> \mathrm{hist})\\ for \\k = 1, \dots, \mathrm{future}\\.

## Usage

``` r
G2G_static_cond_survival(solution, X_vect, hist, future)
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

- hist:

  Number of periods the subject has survived without an event, a whole
  number \>= 0. `hist = 0` gives the unconditional survival curve.

- future:

  Number of periods to forecast ahead, a whole number \>= 1.

## Value

Numeric vector of length `future + 1`: a leading 1 (no additional
periods), then the probability of surviving 1, 2, ..., `future`
additional periods given survival through `hist` periods.

## Details

With the constant covariate effect \\c = \exp(x'\beta)\\, surviving
`hist` periods updates the baseline rate to \\\mathrm{Gamma}(r, \alpha +
\mathrm{hist} \cdot c)\\, and averaging the future survival over that
posterior gives \$\$P(T \> \mathrm{hist} + k \mid T \> \mathrm{hist}) =
\left(\frac{\alpha + \mathrm{hist}\\ c} {\alpha + (\mathrm{hist} + k)\\
c}\right)^r,\$\$ which equals \\S(\mathrm{hist} + k) /
S(\mathrm{hist})\\ from
[`G2G_static_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_survival.md).

## See also

[`G2G_static_posterior_p`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_posterior_p.md),
[`G2G_varying_cond_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_cond_survival.md)
for time-varying covariates.

## Examples

``` r
# \donttest{
  if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
    # Survived 30 days; probability of surviving each of the next 10 days
    G2G_static_cond_survival(fit, data.frame(age = 60, karno = 70),
                             hist = 30, future = 10)
  }
#>  [1] 1.0000000 0.9931628 0.9863809 0.9796537 0.9729807 0.9663614 0.9597953
#>  [8] 0.9532819 0.9468207 0.9404113 0.9340531
# }
```
