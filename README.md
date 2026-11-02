# G2Gcov

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

**G2Gcov** implements the Grassia(II)-Geometric (G2G) distribution for
discrete-time survival analysis with both time-varying and time-invariant
covariates, following the methodology of Fader and Hardie (1997, 2020).

The G2G model accounts for unobserved heterogeneity across subjects through a
Gamma-distributed frailty. Covariates enter through a log-linear link on the
cumulative hazard, making it straightforward to assess whether dynamic factors
— such as marketing interventions or clinical treatment changes — accelerate or
delay events over time.

---

## Installation

Install the released version from CRAN:

```r
install.packages("G2Gcov")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("kaloklee/G2Gcov")
```

---

## Quick Start

### 1. Primary Focus: Time-Varying Covariates (`kb_data_sample`)

For dynamic marketing interventions, subscription tracking, or clinical regimes where covariates evolve over time, data must be structured in **person-period (long) format**. 

Here we demonstrate estimation using `kb_data_sample`, a 20-subject random sample included for rapid experimentation and testing:

```r
library(G2Gcov)

# Load the quick-start sample dataset (20 subjects)
data(kb_data_sample)

# Estimate the discrete-time model with dynamic marketing covariates
fit_varying <- G2G_varying_MLE(
  fo      = "Surv(week, status) ~ coupon + anyp",
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

# Predict period-to-period retention r(t) along the same trajectory
G2G_varying_retention(solution = fit_varying, X_mat = customer_path)
```

> **Note on Datasets:** `kb_data_sample` is a fast-running, 20-subject subset intended for quick trials, unit tests, and interactive exploration. For full empirical analyses, use `data(kb_data)`, which contains the entire cohort of 1,499 subjects (34,506 person-period rows).

---

### 2. Secondary: Time-Invariant (Static) Covariates

When covariates remain fixed throughout observation (such as baseline demographics or treatment arms), standard survival datasets (like `survival::veteran`) can be modeled directly, including both numeric and categorical/factor variables:

```r
library(G2Gcov)
library(survival)

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

# Predict retention over 10 periods
G2G_static_retention(fit_static, data.frame(age = 49, celltype = "squamous"), 10)
```

---

## Data Format

* **Time-invariant models:** One row per subject.
* **Time-varying models:** **Person-period (long) format**, with one row per subject per discrete time period:

| Column | Type | Description |
|---|---|---|
| `id` | Identifier | Unique subject identifier (numeric, character, or factor) |
| `week` (or `time`) | Integer | Sequential discrete period index (1, 2, 3, …) |
| `status` | Binary (0/1) | 1 if the event occurred in that specific period; 0 otherwise |
| `coupon, anyp, ...` | Numeric / Factor / Character | Dynamic or static covariates (continuous, binary, or multi-level categories) |

### Support for Categorical and Non-Numeric Variables

The package natively handles non-numeric variables (both `factor` and `character` vectors) through standard R formula mechanics (`model.matrix`):

* **Automatic Dummy Encoding:** Multi-level factors are automatically expanded into binary dummy indicators using standard treatment contrasts (`contr.treatment`).
* **Baseline Hazard Interpretation:** The first level is treated as the reference category and absorbed into the baseline frailty distribution parameterized by $r$ and $\alpha$ (which represents the hazard when all dummy variables equal zero).
* **Profile Predictions:** When generating predictions with `G2G_static_survival()` or `G2G_static_retention()`, passing a `data.frame` with the original factor levels (e.g., `data.frame(age = 49, celltype = "squamous")`) ensures clean, automatic level mapping.

---

## Main Functions

| Function | Purpose |
|---|---|
| `G2G_varying_MLE()` | Maximum likelihood estimation with time-varying covariates |
| `G2G_static_MLE()` | Maximum likelihood estimation with time-invariant covariates |
| `print()` | S3 print method displaying estimates, SEs, Wald confidence intervals, and AIC |
| `G2G_varying_survival()` | Compute discrete-time survival probabilities $S(t)$ along dynamic covariate paths |
| `G2G_varying_retention()` | Compute period-to-period retention probabilities along dynamic covariate paths |
| `G2G_static_survival()` | Compute discrete-time survival probabilities $S(t)$ under static covariates |
| `G2G_static_retention()` | Compute period-to-period retention probabilities under static covariates |
| `plot_p0()` | Plot baseline per-period event probability density $p_0$ |

---

## Vignettes

To explore full mathematical derivations, likelihood formulations, and empirical comparisons, browse the package vignette:

```r
vignette("G2Gcov", package = "G2Gcov")
```

---

## References

* Fader, P.S. and Hardie, B.G.S. (1997). *How to project customer retention*. **Journal of Interactive Marketing**, 11(1), 61–73.
* Fader, P.S. and Hardie, B.G.S. (2007). *How to project customer retention revisited: The role of duration dependence*. **Journal of Interactive Marketing**, 21(1), 76–90.
* Fader, P.S. and Hardie, B.G.S. (2020). *Incorporating Time-Varying Covariates in a Simple Mixture Model for Discrete-Time Duration Data*. [Research Note](https://brucehardie.com/notes/037/).
* Lee, K.L., Fader, P.S. and Hardie, B.G.S. (2007). *How to project patient persistency*. **Foresight: The International Journal of Applied Forecasting**, Issue 8, 31–35.

---

## License & Contact

Licensed under the [Apache License 2.0](LICENSE).  
Questions, issues, or contributions: Ka Lok Lee (<kaloklee@gmail.com>).