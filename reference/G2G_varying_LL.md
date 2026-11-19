# Log-Likelihood for G2G Model with Time-Varying Covariates

Internal function that computes the negative log-likelihood for the G2G
model with time-varying covariates. Used by the optimization routine.

## Usage

``` r
G2G_varying_LL(par, data_df)
```

## Arguments

- par:

  Numeric vector of parameters on the optimizer's internal scale. The
  first two elements are log(r) and log(alpha) (the logs of the shape
  and rate parameters; both are exponentiated inside the function), and
  the remaining elements are the covariate coefficients.

- data_df:

  Data frame with columns, in this order:

  id

  :   Subject identifier

  time

  :   Time point for each observation. Within each subject, rows must be
      ordered by time, with time running 1, 2, 3, ... and no gaps

  status

  :   Event indicator (0 = no occurrence, 1 = event occurred). A subject
      counts as having had an event if status is 1 on any of their rows,
      and the event is placed in the subject's last period

  ...

  :   Time-varying covariates (remaining columns)

## Value

Scalar value of the negative log-likelihood
