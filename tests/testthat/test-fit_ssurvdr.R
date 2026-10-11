test_that("fit_ssurvdr runs complete workflow with MEC", {
  set.seed(123)

  X <- matrix(rnorm(80 * 10), nrow = 80, ncol = 10)
  time <- rexp(80)
  status <- rbinom(80, 1, 0.7)

  fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 2,
    n_slices = 4,
    covariance = "mec",
    centre = TRUE,
    scale = TRUE
  )

  expect_s3_class(fit, "ssurvdr_fit")
  expect_s3_class(fit$slices, "ssurvdr_slices")
  expect_s3_class(fit$save_fit, "ssurvdr_save_kernel")
  expect_s3_class(fit$basis_fit, "ssurvdr_basis")
  expect_s3_class(fit$cox_fit, "ssurvdr_cox")
  expect_s3_class(fit$icomp_fit, "ssurvdr_icomp")

  expect_true(is.data.frame(fit$criteria))
  expect_equal(ncol(fit$basis_fit$basis), 2)
  expect_true(is.numeric(fit$diagnostics$runtime_seconds))
})

test_that("fit_ssurvdr runs complete workflow with sample covariance", {
  set.seed(123)

  X <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)
  time <- rexp(60)
  status <- rbinom(60, 1, 0.6)

  fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 1,
    n_slices = 3,
    covariance = "sample"
  )

  expect_s3_class(fit, "ssurvdr_fit")
  expect_equal(fit$selected$covariance, "sample")
  expect_equal(fit$selected$d, 1)
})

test_that("fit_ssurvdr rejects invalid d", {
  set.seed(123)

  X <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  time <- rexp(40)
  status <- rbinom(40, 1, 0.5)

  expect_error(
    fit_ssurvdr(
      X = X,
      time = time,
      status = status,
      d = 10
    ),
    "'d' cannot exceed the number of predictors"
  )
})


# test-fit_ssurvdr
set.seed(123)

X <- matrix(rnorm(120 * 15), nrow = 120, ncol = 15)
time <- rexp(120)
status <- rbinom(120, 1, 0.7)

fit <- fit_ssurvdr(
  X = X,
  time = time,
  status = status,
  d = 2,
  n_slices = 4,
  covariance = "mec",
  rho = 0.5,
  alpha = 0.05,
  centre = TRUE,
  scale = TRUE
)
