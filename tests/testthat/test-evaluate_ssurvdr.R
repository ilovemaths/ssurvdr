test_that("evaluate_ssurvdr evaluates fitted model against simulation truth", {
  sim <- simulate_ssurvdr(
    n = 80,
    p = 10,
    active = 1:3,
    covariance = "ar1",
    rho = 0.4,
    seed = 123
  )

  fit <- fit_ssurvdr(
    X = sim$X,
    time = sim$time,
    status = sim$status,
    d = 1,
    n_slices = 3,
    covariance = "mec"
  )

  var_fit <- select_variables(
    basis_fit = fit$basis_fit,
    top = 5,
    variable_names = colnames(sim$X)
  )

  eval_fit <- evaluate_ssurvdr(
    fit = fit,
    simulation = sim,
    variable_selection = var_fit
  )

  expect_s3_class(eval_fit, "ssurvdr_evaluation")
  expect_true(is.data.frame(eval_fit$metrics))
  expect_true(is.data.frame(eval_fit$selected_variables))
  expect_true(is.numeric(eval_fit$projection_distance))
  expect_true(eval_fit$metrics$sensitivity >= 0)
  expect_true(eval_fit$metrics$sensitivity <= 1)
  expect_true(eval_fit$metrics$specificity >= 0)
  expect_true(eval_fit$metrics$specificity <= 1)
})

test_that("evaluate_ssurvdr works with grid object", {
  sim <- simulate_ssurvdr(
    n = 70,
    p = 8,
    active = 1:2,
    covariance = "identity",
    seed = 123
  )

  grid_fit <- fit_ssurvdr_grid(
    X = sim$X,
    time = sim$time,
    status = sim$status,
    d_grid = 1:2,
    covariance_grid = c("mec", "sample"),
    n_slices = 3,
    selector = "ICOMP_IFIM"
  )

  eval_fit <- evaluate_ssurvdr(
    fit = grid_fit,
    simulation = sim,
    top = 4
  )

  expect_s3_class(eval_fit, "ssurvdr_evaluation")
  expect_true(is.data.frame(eval_fit$metrics))
  expect_equal(eval_fit$metrics$n_selected, 4)
})

test_that("evaluate_ssurvdr rejects invalid inputs", {
  sim <- simulate_ssurvdr(
    n = 30,
    p = 5,
    active = 1:2,
    seed = 123
  )

  expect_error(
    evaluate_ssurvdr(
      fit = list(),
      simulation = sim
    ),
    "'fit' must be an ssurvdr_fit or ssurvdr_grid object"
  )

  expect_error(
    evaluate_ssurvdr(
      fit = list(),
      simulation = list()
    ),
    "'simulation' must be an object returned by simulate_ssurvdr"
  )
})

test_that("automatic selection after screening maps back to original predictors", {
  sim <- simulate_ssurvdr(n = 80, p = 10, active = 1:3, seed = 1826)
  retained <- c(7L, 2L, 9L, 1L, 4L)
  fit <- fit_ssurvdr(sim$X[, retained, drop = FALSE], sim$time, sim$status,
                    d = 1, n_slices = 3, covariance = "sample")
  local <- select_variables(fit$basis_fit, top = 2,
                            variable_names = colnames(sim$X)[retained])
  result <- evaluate_ssurvdr(fit, sim, top = 2, screen_indices = retained)
  expect_equal(result$selected_indices, sort(retained[local$selected_indices]))
  expect_true(all(result$selected_indices %in% retained))
})
