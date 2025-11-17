test_that("lppinv works on large underdetermined LP system", {
  skip_if_not_installed("rclsp")
  
  set.seed(123456789)
  
  A_ub <- matrix(rnorm(50 * 500), nrow = 50, ncol = 500)
  A_eq <- matrix(rnorm(25 * 500), nrow = 25, ncol = 500)
  b_ub <- matrix(rnorm(50),       ncol = 1)
  b_eq <- matrix(rnorm(25),       ncol = 1)
  
  result <- lppinv(
    A_ub = A_ub, A_eq = A_eq,
    b_ub = b_ub, b_eq = b_eq,
    non_negative = FALSE,
    final = TRUE,
    alpha = 1.0
  )
  
  # ---- tests ----
  # result is a clsp object, not a wrapper list
  expect_true(is.list(result))
  expect_true("x" %in% names(result))
  
  # x exists and is numeric matrix
  expect_true(is.matrix(result$x))
  
  # stability metrics (if computed)
  expect_true(is.numeric(result$kappaC))
  expect_true(is.numeric(result$kappaB))
  expect_true(is.numeric(result$kappaA))
  
  # NRMSE must be numeric
  expect_true(is.numeric(result$nrmse))
  
  # diagnostic bands must be numeric
  expect_true(all(is.finite(result$x_lower)))
  expect_true(all(is.finite(result$x_upper)))
  
  # bootstrap t-test is available in rclsp
  ttest_res <- rclsp::ttest(
    result,
    sample_size = 30L,
    seed = 123456789L,
    distribution = rnorm,
    partial = TRUE
  )
  
  expect_true(is.list(ttest_res))
  expect_gt(length(ttest_res), 0L)
})
