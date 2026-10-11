test_that("dual sample SAVE reproduces primal SAVE when row space is full", {
  set.seed(1201)

  X <- matrix(rnorm(100 * 8), nrow = 100, ncol = 8)
  time <- rexp(100)
  status <- rbinom(100, 1, 0.7)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 4
  )

  primal <- save_kernel(
    X = X,
    slices = slices,
    covariance = "sample",
    centre = TRUE,
    scale = TRUE,
    backend = "primal"
  )

  dual <- save_kernel(
    X = X,
    slices = slices,
    covariance = "sample",
    centre = TRUE,
    scale = TRUE,
    backend = "dual"
  )

  expect_equal(dual$row_space_rank, 8)
  expect_equal(
    dual$eigenvalues,
    primal$eigenvalues,
    tolerance = 1e-6
  )

  B_primal <- extract_sdr_basis(primal, d = 2)$basis
  B_dual <- extract_sdr_basis(dual, d = 2)$basis

  P_primal <- B_primal %*% t(B_primal)
  P_dual <- B_dual %*% t(B_dual)

  expect_equal(P_dual, P_primal, tolerance = 1e-5)
})


test_that("dual MEC SAVE reproduces primal MEC SAVE when row space is full", {
  set.seed(1202)

  X <- matrix(rnorm(100 * 7), nrow = 100, ncol = 7)
  time <- rexp(100)
  status <- rbinom(100, 1, 0.7)

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = 4
  )

  primal <- save_kernel(
    X = X,
    slices = slices,
    covariance = "mec",
    rho = 0.5,
    alpha = 0.05,
    centre = TRUE,
    scale = TRUE,
    backend = "primal"
  )

  dual <- save_kernel(
    X = X,
    slices = slices,
    covariance = "mec",
    rho = 0.5,
    alpha = 0.05,
    centre = TRUE,
    scale = TRUE,
    backend = "dual"
  )

  expect_equal(dual$row_space_rank, 7)
  expect_equal(
    dual$eigenvalues,
    primal$eigenvalues,
    tolerance = 1e-5
  )

  expect_equal(
    dual$mec_fit$selected_slice_id,
    primal$mec_fit$selected_slice_id
  )

  B_primal <- extract_sdr_basis(primal, d = 2)$basis
  B_dual <- extract_sdr_basis(dual, d = 2)$basis

  P_primal <- B_primal %*% t(B_primal)
  P_dual <- B_dual %*% t(B_dual)

  expect_equal(P_dual, P_primal, tolerance = 1e-4)
})


test_that("dual backend handles p greater than n without dense projection", {
  set.seed(1203)

  X <- matrix(rnorm(60 * 180), nrow = 60, ncol = 180)
  time <- rexp(60)
  status <- rbinom(60, 1, 0.7)

  fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 2,
    n_slices = 3,
    covariance = "mec",
    centre = TRUE,
    scale = TRUE,
    backend = "dual"
  )

  expect_s3_class(fit, "ssurvdr_fit")
  expect_equal(fit$diagnostics$backend, "dual")
  expect_true(fit$diagnostics$row_space_rank < ncol(X))
  expect_equal(dim(fit$basis_fit$basis), c(180, 2))
  expect_null(fit$basis_fit$projection_matrix)
  expect_true(all(is.finite(fit$cox_fit$linear_predictor)))
})


test_that("Cox fitting uses preprocessing learned by SAVE", {
  set.seed(1204)

  X <- matrix(rnorm(90 * 10, mean = 3, sd = 4), nrow = 90, ncol = 10)
  time <- rexp(90)
  status <- rbinom(90, 1, 0.7)

  fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 2,
    n_slices = 4,
    covariance = "sample",
    centre = TRUE,
    scale = TRUE,
    backend = "primal"
  )

  prep <- fit$basis_fit$preprocessing
  X_expected <- sweep(X, 2, prep$center, FUN = "-")
  X_expected <- sweep(X_expected, 2, prep$scale, FUN = "/")
  Z_expected <- X_expected %*% fit$basis_fit$basis
  colnames(Z_expected) <- colnames(fit$cox_fit$reduced_predictors)

  expect_equal(
    fit$cox_fit$reduced_predictors,
    Z_expected,
    tolerance = 1e-10
  )
})


test_that("grid search reuses one SAVE kernel per covariance method", {
  set.seed(1205)

  X <- matrix(rnorm(80 * 12), nrow = 80, ncol = 12)
  time <- rexp(80)
  status <- rbinom(80, 1, 0.7)

  fit <- fit_ssurvdr_grid(
    X = X,
    time = time,
    status = status,
    d_grid = 1:3,
    covariance_grid = c("mec", "sample"),
    n_slices = 3,
    backend = "primal",
    reuse_kernel = TRUE
  )

  expect_true(fit$diagnostics$reuse_kernel)
  expect_equal(length(fit$kernel_cache), 2)
  expect_s3_class(fit$kernel_cache$mec, "ssurvdr_save_kernel")
  expect_s3_class(fit$kernel_cache$sample, "ssurvdr_save_kernel")
  expect_true(all(fit$results$status == "ok"))
})
