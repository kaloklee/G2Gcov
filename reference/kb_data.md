# Kiwi Bubbles Dataset

Person-period discrete-time trial-purchase dataset for "Kiwi Bubbles" (a
masked name for a shelf-stable multipack juice drink aimed primarily at
children) from Market 2 of an IRI BehaviorScan test market, observing
1,499 panelists over 24 weeks alongside weekly market-level
marketing-mix covariates.

## Usage

``` r
kb_data
```

## Format

A data frame with 34,506 rows and 7 variables:

- id:

  Integer panelist identifier.

- week:

  Integer week index for the observation (running 1, 2, ...).

- num_weeks:

  Total weeks observed for the panelist (integer).

- censor:

  Censoring indicator on the panelist's final row: 1 if the panelist was
  censored without a trial purchase, 0 if a trial purchase occurred. 0
  on all preceding rows.

- coupon:

  Exponentially smoothed coupon stock variable capturing carryover
  effects of market-level coupon activity.

- anyp:

  Market-level promotional activity measured as %ACV with any promotion
  (the standard scanner-data sales-weighted promotional measure).

- status:

  Discrete-time trial-purchase event indicator: 1 in the terminal period
  if the panelist made a first purchase, 0 otherwise.

## Source

IRI BehaviorScan "Kiwi Bubbles" test-market dataset (Market 2, first 24
weeks), originally distributed at
<https://brucehardie.com/datasets/kiwibubbles.zip> and processed via
`data-raw/prepare_kb_data.R`.

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
