test_that("save_kernel works with sample covariance", {
  set.seed(123)

  X <- matrix(rnorm(40 * 6), nrow = 40, ncol = 6)

  time <- rexp(40)
  status <- rbinom(40, 1, 0.7)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 4
  )

  fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = "sample"
  )

  expect_s3_class(fit, "ssurvdr_save_kernel")

  expect_true(is.matrix(fit$save_matrix))

  expect_equal(
    dim(fit$save_matrix),
    c(6, 6)
  )

  expect_true(all(is.finite(fit$save_matrix)))
})

test_that("save_kernel works with MEC covariance", {
  set.seed(123)

  X <- matrix(rnorm(50 * 8), nrow = 50, ncol = 8)

  time <- rexp(50)
  status <- rbinom(50, 1, 0.6)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 3
  )

  fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = "mec",
    rho = 0.5,
    alpha = 0.05
  )

  expect_s3_class(fit, "ssurvdr_save_kernel")

  expect_equal(
    dim(fit$save_matrix),
    c(8, 8)
  )

  expect_true("mec_fit" %in% names(fit))
})

test_that("save_kernel rejects invalid covariance option", {
  X <- matrix(rnorm(20), nrow = 10)

  slices <- rep(1:2, each = 5)

  expect_error(
    save_kernel(
      X = X,
      slices = slices,
      covariance = "invalid"
    )
  )
})
