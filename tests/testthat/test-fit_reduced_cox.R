test_that("fit_reduced_cox fits reduced Cox model correctly", {
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

  basis_fit <- extract_sdr_basis(
    save_fit = save_fit,
    d = 2
  )

  cox_fit <- fit_reduced_cox(
    X = X,
    time = time,
    status = status,
    basis_fit = basis_fit
  )

  expect_s3_class(cox_fit, "ssurvdr_cox")

  expect_equal(
    ncol(cox_fit$reduced_predictors),
    2
  )

  expect_true(
    inherits(cox_fit$cox_model, "coxph")
  )

  expect_true(
    is.numeric(cox_fit$aic)
  )

  expect_equal(
    length(cox_fit$coefficients),
    2
  )
})

test_that("fit_reduced_cox rejects invalid status", {
  set.seed(123)

  X <- matrix(rnorm(20 * 4), nrow = 20)

  time <- rexp(20)

  save_fit <- save_kernel(
    X = X,
    slices = rep(1:2, each = 10)
  )

  basis_fit <- extract_sdr_basis(
    save_fit = save_fit,
    d = 2
  )

  expect_error(
    fit_reduced_cox(
      X = X,
      time = time,
      status = rep(2, 20),
      basis_fit = basis_fit
    ),
    "'status' must contain only 0/1 values"
  )
})

test_that("fit_reduced_cox rejects incompatible dimensions", {
  set.seed(123)

  X1 <- matrix(rnorm(40 * 6), nrow = 40)

  X2 <- matrix(rnorm(40 * 5), nrow = 40)

  time <- rexp(40)
  status <- rbinom(40, 1, 0.5)

  save_fit <- save_kernel(
    X = X1,
    slices = rep(1:4, each = 10)
  )

  basis_fit <- extract_sdr_basis(
    save_fit = save_fit,
    d = 2
  )

  expect_error(
    fit_reduced_cox(
      X = X2,
      time = time,
      status = status,
      basis_fit = basis_fit
    ),
    "Number of predictors in X must equal number of rows in basis"
  )
})


# Test-fit_reduced_cox
set.seed(123)

X <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)

time <- rexp(100)
status <- rbinom(100, 1, 0.7)

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

basis_fit <- extract_sdr_basis(
  save_fit = save_fit,
  d = 2
)

cox_fit <- fit_reduced_cox(
  X = X,
  time = time,
  status = status,
  basis_fit = basis_fit
)
