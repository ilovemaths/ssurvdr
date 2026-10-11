test_that("sparsify_sdr_basis preserves exact row support and orthonormality", {
  set.seed(1001)

  p <- 25L
  d <- 3L

  raw_basis <- matrix(
    rnorm(p * d),
    nrow = p,
    ncol = d
  )

  dense_basis <- qr.Q(
    qr(raw_basis),
    complete = FALSE
  )

  basis_fit <- list(
    basis = dense_basis,
    projection_matrix = NULL,
    eigenvalues = seq_len(d),
    full_eigenvalues = seq_len(d),
    explained_ratio = rep(1 / d, d),
    cumulative_ratio = cumsum(rep(1 / d, d)),
    diagnostics = data.frame(component = seq_len(d)),
    d = d,
    p = p,
    covariance_method = "sample",
    preprocessing = list(center = rep(0, p), scale = NULL),
    backend = "dual",
    data_span_restricted = TRUE,
    row_space_rank = 10L,
    null_eigenvalue = 1,
    null_multiplicity = p - 10L,
    n_row_space_eigenvalues_above_null = 3L,
    projection_stored = FALSE
  )

  class(basis_fit) <- "ssurvdr_basis"

  sparse_fit <- sparsify_sdr_basis(
    basis_fit = basis_fit,
    support_size = 8L
  )

  expect_s3_class(sparse_fit, "ssurvdr_sparse_basis")
  expect_true(inherits(sparse_fit, "ssurvdr_basis"))
  expect_equal(sparse_fit$support_size, 8L)
  expect_equal(length(sparse_fit$selected_indices), 8L)

  excluded <- setdiff(
    seq_len(p),
    sparse_fit$selected_indices
  )

  expect_true(
    all(sparse_fit$basis[excluded, , drop = FALSE] == 0)
  )

  expect_equal(
    crossprod(sparse_fit$basis),
    diag(d),
    tolerance = 1e-8
  )

  expected_top <- order(
    sqrt(rowSums(dense_basis^2)),
    decreasing = TRUE
  )[seq_len(8L)]

  expect_setequal(
    sparse_fit$selected_indices,
    expected_top
  )

  expect_true(
    sparse_fit$retained_energy_proportion > 0 &&
      sparse_fit$retained_energy_proportion <= 1
  )

  expect_true(
    sparse_fit$projection_distance_from_dense >= 0
  )
})


test_that("full support reproduces the dense SDR subspace", {
  set.seed(1002)

  p <- 20L
  d <- 2L

  B <- qr.Q(
    qr(matrix(rnorm(p * d), p, d)),
    complete = FALSE
  )

  basis_fit <- list(
    basis = B,
    projection_matrix = NULL,
    eigenvalues = c(2, 1),
    full_eigenvalues = c(2, 1),
    explained_ratio = c(2 / 3, 1 / 3),
    cumulative_ratio = c(2 / 3, 1),
    diagnostics = data.frame(component = 1:2),
    d = d,
    p = p,
    covariance_method = "mec",
    preprocessing = list(center = rep(0, p), scale = NULL),
    backend = "dual",
    data_span_restricted = TRUE,
    row_space_rank = 10L,
    null_eigenvalue = 1,
    null_multiplicity = 10L,
    n_row_space_eigenvalues_above_null = 2L,
    projection_stored = FALSE
  )

  class(basis_fit) <- "ssurvdr_basis"

  full_fit <- sparsify_sdr_basis(
    basis_fit = basis_fit,
    support_size = p
  )

  expect_equal(
    full_fit$retained_energy_proportion,
    1,
    tolerance = 1e-10
  )

  expect_identical(
    full_fit$projection_distance_from_dense,
    0
  )

  expect_identical(
    full_fit$basis,
    B
  )

  expect_true(
    isTRUE(full_fit$full_support_identity)
  )

  expect_equal(
    crossprod(full_fit$basis),
    diag(d),
    tolerance = 1e-8
  )
})


test_that("support smaller than d is rejected", {
  set.seed(1003)

  B <- qr.Q(
    qr(matrix(rnorm(30), 10, 3)),
    complete = FALSE
  )

  basis_fit <- list(
    basis = B,
    d = 3L,
    p = 10L,
    preprocessing = list(center = rep(0, 10), scale = NULL)
  )

  class(basis_fit) <- "ssurvdr_basis"

  expect_error(
    sparsify_sdr_basis(
      basis_fit = basis_fit,
      support_size = 2L
    ),
    "at least the structural dimension"
  )
})
