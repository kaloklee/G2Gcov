# Numerically Stable Log(exp(a) - exp(b))

Computes log(exp(a) - exp(b)) in a numerically stable way, avoiding
overflow and underflow issues. Follows the algorithm of Mächler (2012).

## Usage

``` r
logdiffexp(a, b)
```

## Arguments

- a:

  Numeric vector representing the first value

- b:

  Numeric vector representing the second value

## Value

Numeric vector. Returns -Inf when a == b (including a == b == -Inf), and
NaN when a \< b or a == Inf

## References

Mächler, M. (2012). Accurately Computing log(1 - exp(-\|a\|)).
<https://CRAN.R-project.org/package=Rmpfr>
