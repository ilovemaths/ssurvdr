test_that("simulate_ssurvdr simulates valid identity covariance data", {
  sim <- simulate_ssurvdr(
    n = 50,
    p = 10,
    active = 1:3,
    covariance = "identity",
    seed = 123
  )

  expect_s3_class(sim, "ssurvdr_simulation")
  expect_equal(dim(sim$X), c(50, 10))
  expect_equal(length(sim$time), 50)
  expect_equal(length(sim$status), 50)
  expect_true(all(sim$status %in% c(0L, 1L)))
  expect_equal(dim(sim$true_basis), c(10, 1))
  expect_equal(sim$active, 1:3)
  expect_true(sim$censoring_proportion >= 0)
  expect_true(sim$censoring_proportion <= 1)
})

test_that("simulate_ssurvdr simulates AR1 covariance data", {
  sim <- simulate_ssurvdr(
    n = 40,
    p = 8,
    active = c(1, 4),
    covariance = "ar1",
    rho = 0.4,
    seed = 123
  )

  expect_equal(dim(sim$Sigma), c(8, 8))
  expect_equal(sim$Sigma[1, 2], 0.4)
  expect_equal(sim$Sigma[1, 3], 0.4^2)
})

test_that("simulate_ssurvdr simulates compound covariance data", {
  sim <- simulate_ssurvdr(
    n = 40,
    p = 8,
    active = c(2, 5),
    covariance = "compound",
    rho = 0.3,
    seed = 123
  )

  expect_equal(diag(sim$Sigma), rep(1, 8))
  expect_equal(sim$Sigma[1, 2], 0.3)
})

test_that("simulate_ssurvdr simulates block covariance data", {
  sim <- simulate_ssurvdr(
    n = 40,
    p = 10,
    active = c(1, 2),
    covariance = "block",
    rho = 0.5,
    block_size = 5,
    seed = 123
  )

  expect_equal(sim$Sigma[1, 2], 0.5)
  expect_equal(sim$Sigma[1, 6], 0)
})

test_that("simulate_ssurvdr rejects invalid active indices", {
  expect_error(
    simulate_ssurvdr(
      n = 30,
      p = 5,
      active = c(1, 6)
    ),
    "All entries of 'active' must be valid predictor indices"
  )
})
