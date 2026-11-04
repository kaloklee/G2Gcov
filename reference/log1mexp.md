# Numerically Stable Log(1 - exp(x))

Computes log(1 - exp(x)) in a numerically stable way following the
algorithm of Mächler (2012).

## Usage

``` r
log1mexp(x)
```

## Arguments

- x:

  Numeric vector of values (should be \<= 0 for valid results)

## Value

Numeric vector. Returns -Inf when x == 0 and NaN when x \> 0

## References

Mächler, M. (2012). Accurately Computing log(1 - exp(-\|a\|)).
<https://CRAN.R-project.org/package=Rmpfr>
