test_that("sparse support path contains the dense benchmark and returns a valid fit", {
  skip_if_not_installed("survival")

  set.seed(2001)

  n <- 70L
  p <- 24L

  X <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )

  colnames(X) <- paste0("X", seq_len(p))

  eta <- 0.8 * X[, 1] - 0.5 * X[, 2] + 0.3 * X[, 3]
  event_time <- rexp(n, rate = exp(eta) / 5)
  censor_time <- rexp(n, rate = 0.12)
  time <- pmax(pmin(event_time, censor_time), 1e-4)
  status <- as.integer(event_time <= censor_time)

  dense_fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 2L,
    n_slices = 3L,
    covariance = "sample",
    centre = TRUE,
    scale = TRUE,
    backend = "primal"
  )

  path_fit <- fit_ssurvdr_sparse_path(
    X = X,
    time = time,
    status = status,
    dense_fit = dense_fit,
    support_grid = c(4L, 8L, 12L, p),
    selector = "ICOMP_IFIM"
  )

  expect_s3_class(path_fit, "ssurvdr_sparse_path")
  expect_true(any(path_fit$path$is_dense_support))
  expect_equal(path_fit$d, 2L)
  expect_equal(path_fit$p, p)
  expect_equal(path_fit$diagnostics$n_candidates, 4L)
  expect_true(path_fit$diagnostics$n_successful >= 1L)
  expect_s3_class(path_fit$selected_fit$basis_fit, "ssurvdr_sparse_basis")
  expect_s3_class(path_fit$selected_fit$cox_fit, "ssurvdr_cox")
  expect_s3_class(path_fit$selected_fit$icomp_fit, "ssurvdr_icomp")
})


test_that("full support sparse-path Cox fit is equivalent to dense Cox fit", {
  skip_if_not_installed("survival")

  set.seed(2002)

  n <- 80L
  p <- 18L

  X <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )

  eta <- 0.7 * X[, 1] - 0.4 * X[, 4]
  event_time <- rexp(n, rate = exp(eta) / 6)
  censor_time <- rexp(n, rate = 0.10)
  time <- pmax(pmin(event_time, censor_time), 1e-4)
  status <- as.integer(event_time <= censor_time)

  dense_fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 2L,
    n_slices = 3L,
    covariance = "sample",
    centre = TRUE,
    scale = TRUE,
    backend = "primal"
  )

  path_fit <- fit_ssurvdr_sparse_path(
    X = X,
    time = time,
    status = status,
    dense_fit = dense_fit,
    support_grid = p,
    selector = "ICOMP_IFIM"
  )

  dense_row <- path_fit$dense_row

  expect_equal(
    dense_row$ICOMP_IFIM,
    dense_fit$icomp_fit$icomp_ifim,
    tolerance = 1e-7
  )

  expect_equal(
    dense_row$ICOMP_MISSPEC,
    dense_fit$icomp_fit$icomp_misspec,
    tolerance = 1e-7
  )

  expect_identical(
    path_fit$selected_fit$basis_fit$projection_distance_from_dense,
    0
  )

  expect_true(
    isTRUE(path_fit$selected_fit$basis_fit$full_support_identity)
  )
})


test_that("support path automatically adds the dense support candidate", {
  skip_if_not_installed("survival")

  set.seed(2003)

  n <- 60L
  p <- 15L

  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n) + 0.1
  status <- rbinom(n, 1, 0.65)

  dense_fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 1L,
    n_slices = 3L,
    covariance = "sample",
    backend = "primal"
  )

  path_fit <- fit_ssurvdr_sparse_path(
    X = X,
    time = time,
    status = status,
    dense_fit = dense_fit,
    support_grid = c(3L, 5L, 10L)
  )

  expect_true(p %in% path_fit$support_grid)
  expect_true(any(path_fit$path$is_dense_support))
})
