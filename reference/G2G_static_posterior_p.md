# Posterior Mean of the Per-Period Event Probability Given History (Time-Invariant Covariates)

Updates the population distribution of the baseline rate using the
number of periods a subject has survived, and returns the subject's
posterior mean per-period event probability, which includes the effect
of the subject's own (constant) covariate values.

## Usage

``` r
G2G_static_posterior_p(solution, X_vect, hist)
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
  number \>= 0. `hist = 0` returns the prior mean.

## Value

A single number: the posterior mean of the subject's per-period event
probability given survival through `hist` periods.

## Details

With time-invariant covariates the subject's covariate effect is the
constant \\c = \exp(x'\beta)\\. Surviving `hist` periods gives
\\C(\mathrm{hist}) = \mathrm{hist} \cdot c\\ and the posterior \\\lambda
\sim \mathrm{Gamma}(r, \alpha + \mathrm{hist} \cdot c)\\. The subject's
per-period event probability is \\p = 1 - e^{-\lambda c}\\, and its
posterior mean is \$\$1 - \left(\frac{\alpha + \mathrm{hist}\\ c}
{\alpha + (\mathrm{hist} + 1)\\ c}\right)^r.\$\$ This is also the
probability that the subject has the event in the next period given
survival through `hist` periods, that is, 1 minus the first step of
[`G2G_static_cond_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_cond_survival.md)
with `future = 1`. For a subject whose covariates are all zero (\\c =
1\\) it is the baseline probability \\p_0\\.

## See also

[`G2G_static_cond_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_cond_survival.md),
[`G2G_varying_posterior_p0`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_posterior_p0.md)
for the baseline probability with time-varying covariates.

## Examples

``` r
# \donttest{
  if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    fit <- G2G_static_MLE(Surv(time, status) ~ age + karno, data = veteran)
    # A 60-year-old with Karnofsky score 70 who has survived 30 days
    G2G_static_posterior_p(fit, data.frame(age = 60, karno = 70), hist = 30)
  }
#> [1] 0.006837165
# }
```
