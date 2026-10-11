test_that("compute_icomp computes all criteria", {
  set.seed(123)

  X <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)

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

  icomp_fit <- compute_icomp(
    cox_fit = cox_fit
  )

  expect_s3_class(icomp_fit, "ssurvdr_icomp")

  expect_true(is.data.frame(icomp_fit$criteria))

  expect_true(
    all(c(
      "AIC",
      "BIC",
      "CAIC",
      "ICOMP_IFIM",
      "ICOMP_MISSPEC"
    ) %in% icomp_fit$criteria$criterion)
  )

  expect_true(is.finite(icomp_fit$aic))
  expect_true(is.finite(icomp_fit$bic))
  expect_true(is.finite(icomp_fit$caic))
  expect_true(is.finite(icomp_fit$icomp_ifim))

  expect_true(is.matrix(icomp_fit$vcov_ifim))
  expect_equal(
    dim(icomp_fit$vcov_ifim),
    c(2, 2)
  )
})

test_that("compute_icomp rejects invalid input", {
  expect_error(
    compute_icomp(cox_fit = list()),
    "'cox_fit' must be an object returned by fit_reduced_cox"
  )
})

test_that("compute_icomp rejects invalid eps", {
  set.seed(123)

  X <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)

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
    covariance = "sample"
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

  expect_error(
    compute_icomp(
      cox_fit = cox_fit,
      eps = 0
    ),
    "'eps' must be a single positive numeric value"
  )
})


# test-compute_icomp
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

icomp_fit <- compute_icomp(
  cox_fit = cox_fit
)
