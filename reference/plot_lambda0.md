# Plot the Baseline Rate Distribution of a Fitted G2Gcov Model

Plots the Gamma mixing distribution of the latent baseline rate
\\\lambda_0\\, with all covariates equal to zero, implied by a model
fitted with
[`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md)
(time-invariant covariates) or
[`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md)
(time-varying covariates). The mean and median of \\\lambda_0\\ are
marked with vertical lines and shown in a second title.

## Usage

``` r
plot_lambda0(solution, upper_q = 0.95)
```

## Arguments

- solution:

  A fitted model object of class `"G2Gcov"`.

- upper_q:

  Numeric quantile in `(0, 1)` used to set the upper x-axis limit.
  Default `0.95`.

## Value

Invisibly returns a list with the `mean` and `median` of \\\lambda_0\\.

## Details

In the G2G model, \\p_0 = 1 - e^{-\lambda_0}\\, where \\\lambda_0\\
follows a Gamma distribution with shape \\r\\ and rate \\\alpha\\. The
fitted `r` and `alpha` are on the original covariate scale, so the curve
describes a subject whose covariates are all zero, which can lie outside
the range of the data. The mean and median have closed forms: \\r /
\alpha\\ and `qgamma(0.5, shape = r, rate = alpha)`. When \\r \< 1\\ the
density is unbounded at 0, so the curve rises steeply there and the mean
lies above the median.

## See also

[`plot_p0`](https://kaloklee.github.io/G2Gcov/reference/plot_p0.md),
[`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md),
[`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md)

## Examples

``` r
# \donttest{
  data(kb_data_sample)
  fit <- G2G_varying_MLE(Surv(week, status) ~ coupon + anyp,
                         data = kb_data_sample, subject = "id")
  plot_lambda0(fit)

# }
```
