test_that("select_variables works with ssurvdr_basis object", {
  set.seed(123)

  X <- matrix(rnorm(80 * 10), nrow = 80, ncol = 10)
  colnames(X) <- paste0("gene", seq_len(10))

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

  var_fit <- select_variables(
    basis_fit = basis_fit,
    top = 5,
    variable_names = colnames(X)
  )

  expect_s3_class(var_fit, "ssurvdr_variable_selection")
  expect_true(is.data.frame(var_fit$scores))
  expect_true(is.data.frame(var_fit$selected))
  expect_equal(var_fit$n_selected, 5)
  expect_equal(length(var_fit$selected_variables), 5)
})

test_that("select_variables works with numeric basis matrix", {
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

  var_fit <- select_variables(
    basis_fit = B,
    threshold = 0.25,
    variable_names = paste0("x", 1:4)
  )

  expect_s3_class(var_fit, "ssurvdr_variable_selection")
  expect_equal(var_fit$n_selected, 2)
  expect_true(all(var_fit$selected_variables %in% c("x1", "x4")))
})

test_that("select_variables rejects invalid inputs", {
  B <- matrix(rnorm(10), nrow = 5)

  expect_error(
    select_variables(
      basis_fit = B,
      threshold = -1
    ),
    "'threshold' must be a single non-negative numeric value"
  )

  expect_error(
    select_variables(
      basis_fit = B,
      top = 10
    ),
    "'top' cannot exceed the number of predictors"
  )

  expect_error(
    select_variables(
      basis_fit = B,
      variable_names = c("x1", "x2")
    ),
    "'variable_names' must have length equal to the number of predictors"
  )
})
