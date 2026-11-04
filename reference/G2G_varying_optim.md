# Internal Optimization for Time-Varying Model

Internal wrapper function that performs the actual optimization for the
G2G model with time-varying covariates using the BFGS method. The shape
and rate parameters are optimized on the log scale, so an unconstrained
method is appropriate.

## Usage

``` r
G2G_varying_optim(model_data)
```

## Arguments

- model_data:

  Data frame prepared by `G2G_varying_MLE` containing id, time, status,
  and the (standardized) covariate columns

## Value

The unmodified list returned by
[`optim`](https://rdrr.io/r/stats/optim.html), with elements including:

- par:

  Parameter estimates on the internal scale: log(r), log(alpha) and the
  coefficients of the standardized covariates

- value:

  Negative log-likelihood at the optimum

- convergence:

  Convergence code (0 indicates successful convergence)

- hessian:

  Hessian of the negative log-likelihood at the optimum

Standard errors and confidence intervals are not computed here;
`G2G_varying_MLE` adds them.
