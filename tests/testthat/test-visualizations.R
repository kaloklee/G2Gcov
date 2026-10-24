# Tests for the plotting functions

# --- plot_insample_holdout_paper() on kb_data -------------------------------

out <- plot_insample_holdout_paper(
  fo       = "Surv(week, status) ~ coupon + anyp",
  data     = kb_with_status(),
  subject  = "id",
  split_at = 12L
)

test_that("plot_insample_holdout_paper returns two ggplots and the fit", {
  expect_named(out, c("p_survival", "p_retention", "solution"))
  expect_s3_class(out$p_survival, "ggplot")
  expect_s3_class(out$p_retention, "ggplot")
  expect_s3_class(out$solution, "G2Gcov")
})

test_that("survival curves start at 100% and track Kaplan-Meier on kb_data", {
  surv <- out$p_survival$data
  expect_equal(surv$S_KM[surv$time == 0], 1)
  expect_equal(surv$S_model[surv$time == 0], 1)
  expect_true(all(surv$S_KM >= 0 & surv$S_KM <= 1))
  expect_true(all(surv$S_model >= 0 & surv$S_model <= 1))
  expect_true(all(diff(surv$S_KM) <= 0))
  expect_true(all(diff(surv$S_model) <= 0))
  # 101 of 1,499 kb_data subjects have an event within 24 weeks
  expect_equal(surv$S_KM[surv$time == 24], 1 - 101 / 1499, tolerance = 0.01)
  expect_lt(max(abs(surv$S_KM - surv$S_model)), 0.03)
})

test_that("retention rates are valid probabilities", {
  ret <- out$p_retention$data
  expect_true(all(ret$R_actual >= 0 & ret$R_actual <= 1))
  expect_true(all(ret$R_model >= 0 & ret$R_model <= 1))
})

# --- plot_survival_insample_g2g() / plot_survival_holdout_g2g() -------------

test_that("in-sample and hold-out survivor curves are valid and consistent", {
  panel <- sim_panel()
  p_in  <- plot_survival_insample_g2g("Surv(time, status) ~ x",
                                      train = panel, subject = "id")
  expect_s3_class(p_in, "ggplot")
  s_in <- p_in$data$S_model
  expect_true(all(s_in > 0 & s_in <= 1))
  expect_true(all(diff(s_in) <= 0))

  p_out <- plot_survival_holdout_g2g("Surv(time, status) ~ x",
                                     train = panel[panel$time <= 5, ],
                                     test  = panel[panel$time >  5, ],
                                     subject = "id")
  expect_s3_class(p_out, "ggplot")
  expect_true(all(p_out$data$time > 5))
  expect_true(all(diff(p_out$data$S_model) <= 0))
})

# --- plot_p0() ---------------------------------------------------------------

test_that("plot_p0 draws the distribution of p0 = 1 - exp(-lambda)", {
  fit <- structure(list(par = c(2, 20)), class = c("G2Gcov", "list"))
  p   <- plot_p0(fit)
  expect_s3_class(p, "ggplot")

  d <- p$data
  expect_true(all(d$x >= 0 & d$x < 1))
  # The shaded band runs from the 5% to the 95% quantile of p0, so it holds
  # 90% of the probability.
  band <- d[d$x >= -expm1(-stats::qgamma(0.05, 2, 20)) &
            d$x <= -expm1(-stats::qgamma(0.95, 2, 20)), ]
  area <- sum(diff(band$x) * (head(band$y, -1) + tail(band$y, -1)) / 2)
  expect_equal(area, 0.90, tolerance = 0.01)
})
