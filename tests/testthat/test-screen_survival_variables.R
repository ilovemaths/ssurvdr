test_that("screen_survival_variables screens variables by top rank", {
  sim <- simulate_ssurvdr(
    n = 80,
    p = 12,
    active = 1:3,
    covariance = "ar1",
    signal_strength = 1,
    seed = 123
  )

  screen_fit <- screen_survival_variables(
    X = sim$X,
    time = sim$time,
    status = sim$status,
    top = 5
  )

  expect_s3_class(screen_fit, "ssurvdr_screening")
  expect_true(is.data.frame(screen_fit$scores))
  expect_true(is.data.frame(screen_fit$selected))
  expect_equal(screen_fit$n_selected, 5)
  expect_equal(length(screen_fit$selected_indices), 5)
})

test_that("screen_survival_variables screens variables by p-value", {
  sim <- simulate_ssurvdr(
    n = 80,
    p = 10,
    active = 1:2,
    signal_strength = 1,
    seed = 123
  )

  screen_fit <- screen_survival_variables(
    X = sim$X,
    time = sim$time,
    status = sim$status,
    p_value = 1
  )

  expect_s3_class(screen_fit, "ssurvdr_screening")
  expect_true(screen_fit$n_selected >= 1)
})

test_that("screen_survival_variables rejects invalid inputs", {
  X <- matrix(rnorm(40 * 5), nrow = 40)
  time <- rexp(40)
  status <- rbinom(40, 1, 0.5)

  expect_error(
    screen_survival_variables(
      X = X,
      time = time,
      status = status,
      top = 10
    ),
    "'top' cannot exceed the number of predictors"
  )

  expect_error(
    screen_survival_variables(
      X = X,
      time = time,
      status = rep(2, 40)
    ),
    "'status' must contain only 0/1 values"
  )
})
