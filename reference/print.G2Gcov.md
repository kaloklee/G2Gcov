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
  data(kb_data_sample)
  fit <- G2G_varying_MLE(Surv(week, status) ~ coupon + anyp,
                         data = kb_data_sample, subject = "id")
  print(fit)
#> G2G Model Fit
#> ---------------------------------------------------- 
#> Convergence:       Converged (code 0)
#> Log-likelihood:    -11.1869
#> AIC:               30.3739
#> Parameters:        4
#> 
#> Coefficients:
#>              Estimate Std.Error Lower.95   Upper.95
#> r (shape)      0.0611    0.1047   0.0021     1.7580
#> alpha (rate)   6.3229   31.1821   0.0004 99709.9104
#> coupon        24.9515   22.4348 -19.0199    68.9230
#> anyp          -0.0159    0.0554  -0.1246     0.0927
#> 
#> Note: Wide confidence intervals for r and alpha reflect a flat likelihood ridge (common with small samples, few observed events, or very low unobserved heterogeneity). Interpret r and alpha jointly (e.g., via plot_p0() or plot_lambda0()).
# }
```
