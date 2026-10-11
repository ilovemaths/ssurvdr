test_that("compute_c1f_weights works with basis matrix only", {
  B <- matrix(
    c(
      0.9, 0.1,
      0.2, 0.1,
      0.0, 0.0,
      0.5, 0.5
    ),
    nrow = 4,
    byrow = TRUE
  )

  fit <- compute_c1f_weights(
    basis_fit = B,
    alpha = 1
  )

  expect_s3_class(fit, "ssurvdr_c1f_weights")
  expect_equal(length(fit$weights), 4)
  expect_true(all(is.finite(fit$weights)))
  expect_true(all(fit$weights > 0))
  expect_true(is.data.frame(fit$diagnostics))
})

test_that("compute_c1f_weights works with basis and matrix", {
  set.seed(123)

  B <- matrix(rnorm(10 * 2), nrow = 10, ncol = 2)

  M <- crossprod(matrix(rnorm(50 * 10), nrow = 50, ncol = 10))

  fit <- compute_c1f_weights(
    basis_fit = B,
    matrix = M,
    alpha = 0.5,
    gamma = 1,
    cap = 10
  )

  expect_s3_class(fit, "ssurvdr_c1f_weights")
  expect_equal(length(fit$weights), 10)
  expect_true(is.finite(fit$c1f))
  expect_true(all(fit$weights <= 10))
})

test_that("compute_c1f_weights works with ssurvdr basis object", {
  set.seed(123)

  X <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)
  time <- rexp(60)
  status <- rbinom(60, 1, 0.7)

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
    d = 2
  )

  weight_fit <- compute_c1f_weights(
    basis_fit = basis_fit,
    matrix = save_fit$save_matrix
  )

  expect_s3_class(weight_fit, "ssurvdr_c1f_weights")
  expect_equal(length(weight_fit$weights), 8)
})

test_that("compute_c1f_weights rejects invalid inputs", {
  B <- matrix(rnorm(10), nrow = 5)

  expect_error(
    compute_c1f_weights(
      basis_fit = B,
      alpha = -1
    ),
    "'alpha' must be a single non-negative numeric value"
  )

  expect_error(
    compute_c1f_weights(
      basis_fit = B,
      gamma = 0
    ),
    "'gamma' must be a single positive numeric value"
  )

  expect_error(
    compute_c1f_weights(
      basis_fit = B,
      matrix = diag(4)
    ),
    "'matrix' dimension must match the number of rows in the basis matrix"
  )
})
