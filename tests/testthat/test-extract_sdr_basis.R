test_that("extract_sdr_basis extracts SDR directions correctly", {
  set.seed(123)

  X <- matrix(rnorm(50 * 6), nrow = 50, ncol = 6)

  time <- rexp(50)
  status <- rbinom(50, 1, 0.7)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 4
  )

  save_fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = "sample"
  )

  basis_fit <- extract_sdr_basis(
    save_fit = save_fit,
    d = 2
  )

  expect_s3_class(basis_fit, "ssurvdr_basis")

  expect_equal(
    dim(basis_fit$basis),
    c(6, 2)
  )

  expect_equal(
    dim(basis_fit$projection_matrix),
    c(6, 6)
  )

  expect_true(all(is.finite(basis_fit$basis)))
})

test_that("extract_sdr_basis works with MEC SAVE kernel", {
  set.seed(123)

  X <- matrix(rnorm(40 * 8), nrow = 40, ncol = 8)

  time <- rexp(40)
  status <- rbinom(40, 1, 0.6)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 3
  )

  save_fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = "mec"
  )

  basis_fit <- extract_sdr_basis(
    save_fit = save_fit,
    d = 3
  )

  expect_equal(
    dim(basis_fit$basis),
    c(8, 3)
  )

  expect_true(all(basis_fit$eigenvalues >= 0))
})

test_that("extract_sdr_basis rejects invalid d", {
  set.seed(123)

  X <- matrix(rnorm(20 * 4), nrow = 20, ncol = 4)

  time <- rexp(20)
  status <- rbinom(20, 1, 0.5)

  slices <- survival_slicing(
    time = time,
    status = status
  )

  save_fit <- save_kernel(
    X = X,
    slices = slices
  )

  expect_error(
    extract_sdr_basis(save_fit, d = 10),
    "'d' cannot exceed the number of predictors"
  )
})
