library(survival)

# -------------------------------------------------------------------------
# 1. Dataset Verification
# -------------------------------------------------------------------------

test_that("kb_data_sample dataset is loaded and structured correctly", {
  expect_true(exists("kb_data_sample"))
  expect_true(is.data.frame(kb_data_sample))
  expect_true(all(c("id", "week", "status", "coupon") %in% names(kb_data_sample)))
  expect_gt(nrow(kb_data_sample), 0)
})

# -------------------------------------------------------------------------
# 2. G2G_varying_MLE Tests
# -------------------------------------------------------------------------

test_that("G2G_varying_MLE fits on kb_data_sample and returns correct G2Gcov structure", {
  fit <- G2G_varying_MLE(
    Surv(week, status) ~ coupon,
    data = kb_data_sample,
    subject = "id"
  )
  
  # S3 class inheritance
  expect_s3_class(fit, "G2Gcov")
  
  # Verified top-level components
  expected_slots <- c("par", "value", "counts", "convergence", "hessian",
                      "par_stderr", "par_lower", "par_upper", "center", "scale", "level")
  expect_true(all(expected_slots %in% names(fit)))
  
  # Optimization converged
  expect_equal(fit$convergence, 0)
  
  # Check parameter vector lengths and standard errors
  expect_equal(length(fit$par), length(fit$par_stderr))
  expect_equal(length(fit$par), length(fit$par_lower))
  expect_equal(length(fit$par), length(fit$par_upper))
  
  # Shape parameters r and alpha must be strictly positive (positional 1 and 2)
  expect_gt(unname(fit$par[1]), 0)
  expect_gt(unname(fit$par[2]), 0)
  
  # Parameter estimates and standard errors are finite and non-NA
  expect_false(any(is.na(fit$par)))
  expect_true(all(fit$par_stderr > 0, na.rm = TRUE))
})

test_that("G2G_varying_MLE S3 print method executes cleanly", {
  fit <- G2G_varying_MLE(
    Surv(week, status) ~ coupon,
    data = kb_data_sample,
    subject = "id"
  )
  
  expect_output(print(fit))
})

test_that("G2G_varying_MLE fits multiple covariates (coupon + anyp)", {
  expect_no_error({
    fit_multi <- G2G_varying_MLE(
      Surv(week, status) ~ coupon + anyp,
      data = kb_data_sample,
      subject = "id"
    )
  })
})

test_that("G2G_varying_MLE respects custom Wald confidence level", {
  fit_90 <- G2G_varying_MLE(
    Surv(week, status) ~ coupon,
    data = kb_data_sample,
    subject = "id",
    level = 0.90
  )
  
  fit_99 <- G2G_varying_MLE(
    Surv(week, status) ~ coupon,
    data = kb_data_sample,
    subject = "id",
    level = 0.99
  )
  
  # 99% interval width must be strictly wider than 90% interval width
  width_90 <- fit_90$par_upper - fit_90$par_lower
  width_99 <- fit_99$par_upper - fit_99$par_lower
  expect_true(all(width_99 > width_90))
})

# -------------------------------------------------------------------------
# 3. G2G_static_MLE Tests
# -------------------------------------------------------------------------

test_that("G2G_static_MLE fits and returns correct G2Gcov structure", {
  # Select last observation per customer ID
  kb_static <- do.call(rbind, lapply(split(kb_data_sample, kb_data_sample$id), function(df) {
    df[which.max(df$week), ]
  }))
  
  fit_static <- G2G_static_MLE(
    Surv(week, status) ~ coupon,
    data = kb_static
  )
  
  # S3 class inheritance
  expect_s3_class(fit_static, "G2Gcov")
  
  # Structure and convergence
  expect_true(all(c("par", "value", "convergence", "par_stderr", "par_lower", "par_upper") %in% names(fit_static)))
  expect_equal(fit_static$convergence, 0)
  expect_equal(length(fit_static$par), length(fit_static$par_stderr))
  
  # Shape parameters r and alpha must be strictly positive (positional 1 and 2)
  expect_gt(unname(fit_static$par[1]), 0)
  expect_gt(unname(fit_static$par[2]), 0)
  
  expect_false(any(is.na(fit_static$par)))
  expect_true(all(fit_static$par_stderr > 0, na.rm = TRUE))
})

test_that("G2G_static_MLE S3 print method executes cleanly", {
  kb_static <- do.call(rbind, lapply(split(kb_data_sample, kb_data_sample$id), function(df) {
    df[which.max(df$week), ]
  }))
  
  fit_static <- G2G_static_MLE(
    Surv(week, status) ~ coupon,
    data = kb_static
  )
  
  expect_output(print(fit_static))
})

test_that("G2G_static_MLE stops on invalid discrete time values", {
  kb_static <- do.call(rbind, lapply(split(kb_data_sample, kb_data_sample$id), function(df) {
    df[which.max(df$week), ]
  }))
  
  # Discrete-time duration models require time >= 1
  kb_static$week[1] <- 0
  
  expect_error(
    G2G_static_MLE(
      Surv(week, status) ~ coupon,
      data = kb_static
    )
  )
})