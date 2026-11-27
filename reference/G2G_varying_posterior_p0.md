# Posterior Mean of the Baseline Per-Period Event Probability Given History

Updates the population distribution of the baseline rate using a
customer's observed history and returns the posterior mean of the
baseline per-period event probability \\p_0 = 1 - e^{-\lambda}\\, i.e.
the per-period event probability at covariate value zero (the reference
levels of any factors).

## Usage

``` r
G2G_varying_posterior_p0(solution, X_hist = NULL)
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
  per period (t = 1, ..., T), in any of the forms accepted by `X_mat` in
  [`G2G_varying_survival`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_survival.md).
  The history must contain only event-free periods. `NULL` (no history)
  returns the prior mean.

## Value

A single number: the posterior mean of the baseline per-period event
probability given survival through the full history.

## Details

A customer who survived `t` periods with covariate path \\z(1), \dots,
z(t)\\ has posterior \\\lambda \sim \mathrm{Gamma}(r, \alpha + C(t))\\,
where \\C(t) = \sum\_{j=1}^{t} e^{z(j)\beta'}\\. The posterior mean of
\\p_0\\ is then \$\$1 - \left(\frac{\alpha + C(t)}{\alpha + C(t) +
1}\right)^r.\$\$ The history enters only through \\C(t)\\.

## Examples

``` r
# \donttest{
  data(kb_data_sample)
  fit <- G2G_varying_MLE(
    Surv(week, status) ~ coupon + anyp,
    data    = kb_data_sample,
    subject = "id"
  )
  # A customer who has survived the first 8 weeks
  hist <- kb_data_sample[kb_data_sample$id == 165, c("coupon", "anyp")][1:8, ]
  G2G_varying_posterior_p0(fit, hist)
#> [1] 0.003069842
# }
```
