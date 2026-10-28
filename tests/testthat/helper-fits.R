# Shared fixtures. Fitting kb_data takes close to a minute, so it is fitted
# once per test run and reused.

kb_with_status <- function() {
  kb <- kb_data
  kb$status <- as.integer(kb$week == kb$num_weeks & kb$censor == 0)
  kb
}

kb_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit))
      fit <<- G2G_varying_MLE("Surv(week, status) ~ coupon + anyp",
                              data = kb_with_status(), subject = "id")
    fit
  }
})

# Small simulated person-period panel with one market-level covariate; fits
# in a few seconds.
sim_panel <- function(n = 300, horizon = 8, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(horizon)
  rows <- lapply(seq_len(n), function(i) {
    lambda <- stats::rgamma(1, shape = 2, rate = 20)
    last <- horizon
    event <- 0L
    for (t in seq_len(horizon)) {
      if (stats::runif(1) < 1 - exp(-lambda * exp(0.5 * x[t]))) {
        last <- t
        event <- 1L
        break
      }
    }
    data.frame(id = i, time = seq_len(last), x = x[seq_len(last)],
               status = c(rep(0L, last - 1L), event))
  })
  do.call(rbind, rows)
}
