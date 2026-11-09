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
#> Warning: The Hessian at the solution is not positive definite (or could not be inverted), so standard errors and confidence intervals are NA.
  print(fit)
#> G2G Model Fit
#> ---------------------------------------------------- 
#> Convergence:       Converged (code 0)
#> Log-likelihood:    -5612.1808
#> AIC:               11232.3616
#> Parameters:        4
#> 
#> Coefficients:
#>                Estimate Std.Error Lower.95 Upper.95
#> r (shape)     2996.6875        NA       NA       NA
#> alpha (rate) 40722.0679        NA       NA       NA
#> coupon         -78.2544        NA       NA       NA
#> anyp             0.0020        NA       NA       NA
#> 
#> Note: some standard errors are NA (the Hessian at the optimum is singular or not positive definite).
# }
```
