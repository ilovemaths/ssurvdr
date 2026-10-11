test_that("mec_covariance returns a valid covariance object", {
  set.seed(123)

  X <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  time <- rexp(20)
  status <- rbinom(20, size = 1, prob = 0.6)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 3,
    min_slice_size = 2
  )

  fit <- mec_covariance(
    X = X,
    slices = slices,
    rho = 0.5,
    alpha = 0.05
  )

  expect_s3_class(fit, "ssurvdr_mec_covariance")
  expect_true(is.matrix(fit$covariance))
  expect_equal(nrow(fit$covariance), ncol(X))
  expect_equal(ncol(fit$covariance), ncol(X))
  expect_true(all(is.finite(fit$covariance)))
  expect_true(is.data.frame(fit$entropy_scores))
  expect_true(fit$selected_slice_id %in% fit$entropy_scores$slice_id)
})

test_that("mec_covariance works with direct slice IDs", {
  set.seed(123)

  X <- matrix(rnorm(12 * 4), nrow = 12, ncol = 4)
  slice_id <- rep(1:3, each = 4)

  fit <- mec_covariance(
    X = X,
    slices = slice_id,
    rho = 0.6,
    alpha = 0.1
  )

  expect_s3_class(fit, "ssurvdr_mec_covariance")
  expect_equal(dim(fit$covariance), c(4, 4))
})

test_that("mec_covariance rejects invalid inputs", {
  X <- matrix(rnorm(10), nrow = 5, ncol = 2)

  expect_error(
    mec_covariance(X = X, slices = c(1, 1, 2)),
    "'slices' must have length equal to nrow\\(X\\)"
  )

  expect_error(
    mec_covariance(X = X, slices = rep(1, 5), rho = 2),
    "'rho' must be a single numeric value in \\[0, 1\\]"
  )

  expect_error(
    mec_covariance(X = X, slices = rep(1, 5), alpha = -0.1),
    "'alpha' must be a single numeric value in \\[0, 1\\]"
  )
})
