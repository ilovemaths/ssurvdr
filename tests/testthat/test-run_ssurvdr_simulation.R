test_that("run_ssurvdr_simulation runs small Monte Carlo experiment", {
  sim_run <- run_ssurvdr_simulation(
    n = 40,
    p = 10,
    active = 1:2,
    covariance = "ar1",
    rho = 0.3,
    signal_strength = 1,
    censoring_rate = 0.1,
    n_rep = 2,
    screening = TRUE,
    screen_top = 6,
    covariance_methods = c("mec", "sample"),
    d_grid = 1:2,
    n_slices = 3,
    variable_top = 4,
    seed = 123
  )

  expect_s3_class(sim_run, "ssurvdr_simulation_run")
  expect_true(is.data.frame(sim_run$results))
  expect_true(is.data.frame(sim_run$summary))
  expect_true(is.data.frame(sim_run$method_selection))
  expect_equal(nrow(sim_run$results), 2)
  expect_equal(sim_run$settings$dimension_regime, "low-dimensional")
})

test_that("run_ssurvdr_simulation labels ultra-high-dimensional regime", {
  sim_run <- run_ssurvdr_simulation(
    n = 20,
    p = 50,
    active = 1:2,
    covariance = "identity",
    n_rep = 1,
    screening = TRUE,
    screen_top = 8,
    covariance_methods = "mec",
    d_grid = 1,
    n_slices = 2,
    variable_top = 4,
    seed = 123
  )

  expect_equal(sim_run$settings$dimension_regime, "ultra-high-dimensional")
})
