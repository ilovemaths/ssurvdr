test_that("summarise_ssurvdr summarises ssurvdr_fit object", {
  set.seed(123)

  X <- matrix(rnorm(80 * 8), nrow = 80, ncol = 8)
  time <- rexp(80)
  status <- rbinom(80, 1, 0.7)

  fit <- fit_ssurvdr(
    X = X,
    time = time,
    status = status,
    d = 2,
    n_slices = 3,
    covariance = "mec"
  )

  summary_fit <- summarise_ssurvdr(fit)

  expect_s3_class(summary_fit, "ssurvdr_summary")
  expect_true(is.data.frame(summary_fit$model_summary))
  expect_true(is.data.frame(summary_fit$criteria))
  expect_true(is.data.frame(summary_fit$slice_table))
})

test_that("summarise_ssurvdr summarises ssurvdr_grid object", {
  set.seed(123)

  X <- matrix(rnorm(70 * 6), nrow = 70, ncol = 6)
  time <- rexp(70)
  status <- rbinom(70, 1, 0.6)

  grid_fit <- fit_ssurvdr_grid(
    X = X,
    time = time,
    status = status,
    d_grid = 1:2,
    covariance_grid = c("mec", "sample"),
    n_slices = 3,
    selector = "ICOMP_IFIM"
  )

  summary_fit <- summarise_ssurvdr(grid_fit)

  expect_s3_class(summary_fit, "ssurvdr_summary")
  expect_true(is.data.frame(summary_fit$best_summary))
  expect_true(is.data.frame(summary_fit$ranked_results))
  expect_true(is.data.frame(summary_fit$all_results))
})

test_that("summarise_ssurvdr summarises variable selection object", {
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

  summary_fit <- summarise_ssurvdr(
    object = var_fit,
    top = 3
  )

  expect_s3_class(summary_fit, "ssurvdr_summary")
  expect_true(is.data.frame(summary_fit$selection_summary))
  expect_equal(nrow(summary_fit$top_variables), 3)
})

test_that("summarise_ssurvdr rejects invalid object", {
  expect_error(
    summarise_ssurvdr(list()),
    "'object' must be an ssurvdr_fit, ssurvdr_grid, or ssurvdr_variable_selection object"
  )
})
