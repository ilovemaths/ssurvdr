test_that("fit_ssurvdr_grid runs model comparison", {
  set.seed(123)

  X <- matrix(rnorm(80 * 8), nrow = 80, ncol = 8)
  time <- rexp(80)
  status <- rbinom(80, 1, 0.7)

  grid_fit <- fit_ssurvdr_grid(
    X = X,
    time = time,
    status = status,
    d_grid = 1:2,
    covariance_grid = c("mec", "sample"),
    n_slices = 3,
    selector = "ICOMP_IFIM",
    centre = TRUE,
    scale = TRUE
  )

  expect_s3_class(grid_fit, "ssurvdr_grid")
  expect_true(is.data.frame(grid_fit$results))
  expect_true(is.data.frame(grid_fit$ranked_results))
  expect_s3_class(grid_fit$best_model, "ssurvdr_fit")
  expect_equal(grid_fit$diagnostics$n_models, 4)
  expect_true(grid_fit$diagnostics$n_successful >= 1)
  expect_equal(grid_fit$selector, "ICOMP_IFIM")
})

test_that("fit_ssurvdr_grid works with one covariance method", {
  set.seed(123)

  X <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)
  time <- rexp(60)
  status <- rbinom(60, 1, 0.6)

  grid_fit <- fit_ssurvdr_grid(
    X = X,
    time = time,
    status = status,
    d_grid = 1:2,
    covariance_grid = "sample",
    n_slices = 3,
    selector = "BIC"
  )

  expect_s3_class(grid_fit, "ssurvdr_grid")
  expect_equal(nrow(grid_fit$results), 2)
  expect_equal(grid_fit$selector, "BIC")
})

test_that("fit_ssurvdr_grid rejects invalid covariance", {
  set.seed(123)

  X <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
  time <- rexp(40)
  status <- rbinom(40, 1, 0.7)

  expect_error(
    fit_ssurvdr_grid(
      X = X,
      time = time,
      status = status,
      covariance_grid = "bad"
    ),
    "'covariance_grid' currently supports only 'mec' and 'sample'"
  )
})
