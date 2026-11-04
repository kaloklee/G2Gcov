# KB Dataset

Person-period discrete-time survival dataset observing 1,499 subjects
for up to 24 weeks with marketing covariates.

## Usage

``` r
kb_data
```

## Format

A data frame with 34,506 rows and 7 variables:

- id:

  Integer subject identifier.

- week:

  Integer week index for the observation (running 1, 2, ...).

- num_weeks:

  Total weeks observed for the subject (integer).

- censor:

  Censoring indicator on the subject's final row: 1 if the subject was
  censored before experiencing the event, 0 if the event occurred. 0 on
  all preceding rows.

- coupon:

  Coupon stock covariate measuring aggregate distributed coupons.

- anyp:

  Trade promotions covariate measuring percentage of stores with active
  promotions in the market during the week.

- status:

  Discrete-time failure event indicator: 1 in the terminal period if the
  subject experienced the event, 0 otherwise.

## Source

Simulated customer churn data modeled after retail cohort dynamics.
Processed via `data-raw/prepare_kb_data.R`.

## Examples

``` r
data(kb_data)
head(kb_data)
#>   id week num_weeks censor coupon  anyp status
#> 1  1    1         1      0      0 55.61      1
#> 2  2    1         1      0      0 55.61      1
#> 3  3    1         1      0      0 55.61      1
#> 4  4    1         1      0      0 55.61      1
#> 5  5    1         1      0      0 55.61      1
#> 6  6    1         1      0      0 55.61      1
```
