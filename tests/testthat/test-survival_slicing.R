test_that("survival_slicing creates valid slices", {
  time <- c(5, 8, 3, 10, 12, 4, 6, 9, 2, 7)
  status <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0)

  fit <- survival_slicing(
    time = time,
    status = status,
    n_slices = 2,
    min_slice_size = 2
  )

  expect_s3_class(fit, "ssurvdr_slices")
  expect_equal(length(fit$slice_id), length(time))
  expect_equal(length(fit$slice_label), length(time))
  expect_false(anyNA(fit$slice_id))
  expect_true(is.data.frame(fit$slice_table))
  expect_true(all(fit$status %in% c(0L, 1L)))
})

test_that("survival_slicing rejects invalid inputs", {
  expect_error(
    survival_slicing(time = c(1, 2, NA), status = c(1, 0, 1)),
    "'time' must not contain missing values"
  )

  expect_error(
    survival_slicing(time = c(1, 2, 3), status = c(1, 0)),
    "'time' and 'status' must have the same length"
  )

  expect_error(
    survival_slicing(time = c(1, 2, 3), status = c(1, 2, 0)),
    "'status' must contain only 0/1 values"
  )
})
