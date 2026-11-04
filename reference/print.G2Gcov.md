# Print a Fitted G2Gcov Model

Displays a formatted summary of a G2G model fit returned by
[`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md)
(time-invariant covariates) or
[`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md)
(time-varying covariates), including convergence status, log-likelihood,
AIC, and a parameter table with standard errors and Wald confidence
intervals at the fitted confidence level.

## Usage

``` r
# S3 method for class 'G2Gcov'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A fitted model object of class `"G2Gcov"`.

- digits:

  Integer. Number of decimal places to print. Default 4.

- ...:

  Further arguments passed to or from other methods (unused).

## Value

Invisibly returns `x`.

## Details

Both fitting functions return an object of class `"G2Gcov"`, and this
one method prints either. It only relies on the elements the two have in
common: `par`, `par_stderr`, `par_lower`, `par_upper`, `value`,
`convergence`, and optional `level`. If any of the standard errors or
confidence limits are missing or `NA` (the fitting functions return `NA`
when the Hessian at the optimum is singular or not positive definite),
they print as `NA` and a note is added below the table.

## See also

[`G2G_static_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_static_MLE.md),
[`G2G_varying_MLE`](https://kaloklee.github.io/G2Gcov/reference/G2G_varying_MLE.md)

## Examples

``` r
# \donttest{
  data(kb_data)
  kb_data$status <- 1L - kb_data$censor
  fit <- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                          data = kb_data, subject = "id")
  print(fit)
#> G2G Model Fit
#> ---------------------------------------------------- 
#> Convergence:       Converged (code 0)
#> Log-likelihood:    -5612.181
#> AIC:               11232.3619
#> Parameters:        4
#> 
#> Coefficients:
#>                Estimate  Std.Error  Lower.95    Upper.95
#> r (shape)     2994.7624  3725.7011  261.4605  34301.9392
#> alpha (rate) 40695.9287 50674.4945 3545.1645 467159.8788
#> coupon         -78.2543     4.2463  -86.5769    -69.9318
#> anyp             0.0020     0.0012   -0.0002      0.0043
# }
```
