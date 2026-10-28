# G2Gcov 0.1.0

* Initial CRAN release.
* `G2G_static_MLE()`: maximum likelihood estimation for G2G models with
  time-invariant covariates. Includes input validation with informative error
  messages for missing columns, NA values, non-binary status, and degenerate
  data. `Surv()` in the formula is read, not evaluated, so the survival
  package does not need to be attached.
* `G2G_varying_MLE()`: maximum likelihood estimation for G2G models with
  time-varying covariates in person-period format. Includes the same input
  validation as the static version plus checks for binary status values,
  positive time periods, the presence of at least one event, and events
  occurring only in a subject's final period. Rows may be supplied in any
  order.
* Both estimators report `r` and `alpha` with delta-method standard errors on
  the natural scale and 95% intervals computed on the log scale.
* `print.G2Gcov()`: S3 print method for fitted model objects. Displays
  convergence status, log-likelihood, AIC, and a formatted parameter table
  with standard errors and 95% Wald confidence intervals.
* `plot_p0()`: visualize the distribution of the baseline per-period event
  probability `p0 = 1 - exp(-lambda)` implied by the fitted Gamma
  distribution of `lambda`.
* `plot_survival_insample_g2g()`: in-sample G2G survivor curve.
* `plot_survival_holdout_g2g()`: hold-out G2G survivor curve fitted on
  training data.
* `plot_insample_holdout_paper()`: paper-style actual-vs-model survival and
  retention plots with optional group faceting.
* The model survivor curves chain together the mean model retention of the
  subjects at risk in each period, mirroring the Kaplan-Meier estimator.
* Includes `kb_data`, a sample person-period dataset used in the package
  vignettes, and `kb_covars_sample`, a small subset of it.
