# G2Gcov

[![License: Apache
2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

**G2Gcov** implements the Grassia(II)-Geometric (G2G) distribution for
discrete-time survival analysis with both time-varying and
time-invariant covariates, following the methodology of Fader et
al. (2019) and Fader and Hardie (2020).

## Installation

Install the released version from CRAN:

``` r

install.packages("G2Gcov")
```

Or install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("kaloklee/G2Gcov")
```

------------------------------------------------------------------------

## Quick Start

### 1. Primary Focus: Time-Varying Covariates (`kb_data_sample`)

For dynamic marketing interventions, subscription tracking, or clinical
regimes where covariates evolve over time, data must be structured in
**person-period (long) format**.

Here we demonstrate estimation using `kb_data_sample`, a 20-subject
random sample included for rapid experimentation and testing:

``` r

library(G2Gcov)

# Load the quick-start sample dataset (20 subjects)
kb_data_sample<-G2Gcov::kb_data_sample

# Estimate the discrete-time model with dynamic marketing covariates
fit_varying <- G2G_varying_MLE(
  Surv(week, status) ~ coupon + anyp,
  data    = kb_data_sample,
  subject = "id"
)

# Print parameter estimates, standard errors, and Wald confidence intervals
print(fit_varying)

# Plot baseline churn probability distribution across the population (p0)
plot_p0(fit_varying)

# Predict individual survival curve S(t) along a specific customer trajectory
# (Subject 24 is included in kb_data_sample)
customer_path <- kb_data_sample[kb_data_sample$id == 24, c("coupon", "anyp")]
G2G_varying_survival(solution = fit_varying, X_mat = customer_path)
```

> **Note on Datasets:** `kb_data_sample` is a fast-running, 20-subject
> subset intended for quick trials, unit tests, and interactive
> exploration. For full empirical analyses, use `data(kb_data)`, which
> contains the entire cohort of 1,499 subjects (34,506 person-period
> rows).

------------------------------------------------------------------------

### 2. Secondary: Time-Invariant (Static) Covariates

When covariates remain fixed throughout observation (such as baseline
demographics or treatment arms), standard survival datasets (like
[`survival::veteran`](https://rdrr.io/pkg/survival/man/veteran.html))
can be modeled directly, including both numeric and categorical/factor
variables:

``` r

library(G2Gcov)
library(survival)

# Load the Veteran dataset
veteran<-survival::veteran

# Estimate the static model with numeric (age) and factor (celltype) covariates
fit_static <- G2G_static_MLE(
  Surv(time, status) ~ age + celltype,
  data = veteran
)

# Inspect model summary
print(fit_static)

# Plot baseline probability distribution (at reference covariate levels)
plot_p0(fit_static)

# Predict survival over 10 periods for a specific baseline covariate profile
G2G_static_survival(fit_static, data.frame(age = 49, celltype = "squamous"), 10)
```

------------------------------------------------------------------------

## Vignettes

For a more comprehensive demonstration of the package’s
functionality,browse the package vignette:

``` r

vignette("introduction_to_g2gcov", package = "G2Gcov")
```

------------------------------------------------------------------------

## References

- Fader, P. S., Hardie, B. G. S., McCarthy, D., & Vaidyanathan, R.
  (2019). Exploring the equivalence of two common mixture models for
  duration data. *The American Statistician*, 73(3), 288–295.
- Fader, P. S. & Hardie, B. G. S. (2020). *Incorporating Time-Varying
  Covariates in a Simple Mixture Model for Discrete-Time Duration Data*.
  [Research Note](https://brucehardie.com/notes/037/).

------------------------------------------------------------------------

## License & Contact

Licensed under the [Apache License
2.0](https://kaloklee.github.io/G2Gcov/LICENSE).  
Questions, issues, or contributions: Ka Lok Lee (<kaloklee@gmail.com>).
