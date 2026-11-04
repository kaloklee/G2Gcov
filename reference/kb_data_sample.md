# KB Sample Dataset

A 20-subject random sample of
[`kb_data`](https://kaloklee.github.io/G2Gcov/reference/kb_data.md)
provided for fast unit testing, minimal package checks, and vignette
examples.

## Usage

``` r
kb_data_sample
```

## Format

A data frame with the same 7 variables as
[`kb_data`](https://kaloklee.github.io/G2Gcov/reference/kb_data.md):

- id:

  Integer subject identifier.

- week:

  Integer week index for the observation.

- num_weeks:

  Total weeks observed for the subject (integer).

- censor:

  Censoring indicator (1 if censored on final row, 0 otherwise).

- coupon:

  Coupon stock covariate.

- anyp:

  Trade promotions covariate.

- status:

  Discrete-time failure event indicator (1 if event occurred in the
  period, 0 otherwise).

## Source

Derived as a random sample from
[`kb_data`](https://kaloklee.github.io/G2Gcov/reference/kb_data.md) via
`data-raw/prepare_kb_data.R`.

## Examples

``` r
data(kb_data_sample)
head(kb_data_sample)
#>   id week num_weeks censor     coupon  anyp status
#> 1 24    1         4      0 0.00000000 55.61      0
#> 2 24    2         4      0 0.00000000 55.82      0
#> 3 24    3         4      0 0.01803704 55.65      0
#> 4 24    4         4      0 0.12625929 81.17      1
#> 5 49    1         7      0 0.00000000 55.61      0
#> 6 49    2         7      0 0.00000000 55.82      0
```
