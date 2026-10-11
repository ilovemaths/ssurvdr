test_that("sparse_projection returns sparse projection object", {
  set.seed(123)

  X <- matrix(rnorm(80 * 10), nrow = 80, ncol = 10)
  time <- rexp(80)
  status <- rbinom(80, 1, 0.7)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 4
  )

  save_fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = "mec"
  )

  sparse_fit <- sparse_projection(
    save_fit = save_fit,
    d = 2,
    lambda = 0.01
  )

  expect_s3_class(sparse_fit, "ssurvdr_sparse_projection")
  expect_equal(dim(sparse_fit$basis), c(10, 2))
  expect_equal(dim(sparse_fit$projection_matrix), c(10, 10))
  expect_equal(dim(sparse_fit$sparse_kernel), c(10, 10))
  expect_true(is.data.frame(sparse_fit$variable_scores))
  expect_true(all(is.finite(sparse_fit$basis)))
})

test_that("sparse_projection works with weights", {
  set.seed(123)

  X <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)
  time <- rexp(60)
  status <- rbinom(60, 1, 0.6)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 3
  )

  save_fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = "sample"
  )

  weights <- seq(1, 2, length.out = 8)

  sparse_fit <- sparse_projection(
    save_fit = save_fit,
    d = 2,
    lambda = 0.02,
    weights = weights
  )

  expect_s3_class(sparse_fit, "ssurvdr_sparse_projection")
  expect_equal(length(sparse_fit$weights), 8)
  expect_equal(dim(sparse_fit$basis), c(8, 2))
})

test_that("sparse_projection rejects invalid inputs", {
  set.seed(123)

  X <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
  time <- rexp(40)
  status <- rbinom(40, 1, 0.7)

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

  expect_error(
    sparse_projection(
      save_fit = save_fit,
      d = 10
    ),
    "'d' cannot exceed the number of predictors"
  )

  expect_error(
    sparse_projection(
      save_fit = save_fit,
      lambda = -1
    ),
    "'lambda' must be a single non-negative numeric value"
  )

  expect_error(
    sparse_projection(
      save_fit = save_fit,
      weights = c(1, 2)
    ),
    "'weights' must be NULL or a numeric vector of length p"
  )
})
