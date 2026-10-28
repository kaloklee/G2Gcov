#' Numerically Stable Log(1 - exp(x))
#'
#' Computes log(1 - exp(x)) in a numerically stable way following the 
#' algorithm of Mächler (2012).
#'
#' @param x Numeric vector of values (should be <= 0 for valid results)
#' @return Numeric vector. Returns -Inf when x == 0 and NaN when x > 0
#' @keywords internal
#' @references 
#' Mächler, M. (2012). Accurately Computing log(1 - exp(-|a|)). 
#' \url{https://CRAN.R-project.org/package=Rmpfr}
log1mexp <- function(x) {
  ifelse(
    x > -0.6931472,
    log(-expm1(x)),
    log1p(-exp(x))
  )
}


#' Numerically Stable Log(exp(a) - exp(b))
#'
#' Computes log(exp(a) - exp(b)) in a numerically stable way, avoiding 
#' overflow and underflow issues. Follows the algorithm of Mächler (2012).
#'
#' @param a Numeric vector representing the first value
#' @param b Numeric vector representing the second value
#' @return Numeric vector. Returns -Inf when a == b (including a == b == -Inf),
#'   and NaN when a < b or a == Inf
#' @keywords internal
#' @references 
#' Mächler, M. (2012). Accurately Computing log(1 - exp(-|a|)). 
#' \url{https://CRAN.R-project.org/package=Rmpfr}
logdiffexp <- function(a, b) {
  ifelse(
    (a < Inf) & (a > b),
    a + log1mexp(b - a),
    ifelse((a < Inf) & (a == b), -Inf, NaN)
  )
}


# r and alpha are optimized as log(r) and log(alpha). Build 95% Wald
# intervals on that scale, then report r and alpha on their natural scale,
# with standard errors from the delta method: SE(exp(x)) = exp(x) * SE(x).
.natural_scale <- function(solution) {
  se  <- solution$par_stderr
  idx <- 1:2
  solution$par_upper <- solution$par + 1.96 * se
  solution$par_lower <- solution$par - 1.96 * se
  solution$par[idx]        <- exp(solution$par[idx])
  solution$par_upper[idx]  <- exp(solution$par_upper[idx])
  solution$par_lower[idx]  <- exp(solution$par_lower[idx])
  solution$par_stderr[idx] <- solution$par[idx] * se[idx]
  solution
}
