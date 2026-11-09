# Check a Number of Periods

Internal helper that stops unless `x` is a single whole number at or
above `min`.

## Usage

``` r
G2G_check_periods(x, name, min)
```

## Arguments

- x:

  Value to check.

- name:

  Argument name used in the error message.

- min:

  Smallest allowed value.

## Value

`NULL`, invisibly, when the value is valid.
