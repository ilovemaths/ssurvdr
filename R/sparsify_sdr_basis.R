#' Row-sparsify an SDR basis
#'
#' Constructs a row-sparse sufficient dimension reduction basis by retaining the
#' predictors with the largest row norms from an existing SDR basis. The retained
#' columns are re-orthonormalised without altering the exact zero pattern of the
#' excluded predictor rows.
#'
#' This function is deliberately a post-estimation row-sparsification operator.
#' It does not claim to solve a penalised SAVE optimisation problem.
#'
#' @param basis_fit Object returned by [extract_sdr_basis()] or another object
#'   inheriting from `ssurvdr_basis`.
#' @param support_size Positive integer number of predictor rows to retain.
#'   Must be at least the structural dimension and no larger than p.
#' @param variable_names Optional character vector of predictor names. If omitted,
#'   row names of the basis are used, followed by `V1`, ..., `Vp` if needed.
#' @param rank_tol Positive numerical tolerance used to verify that the retained
#'   rows still span the requested structural dimension.
#'
#' @return An object inheriting from `ssurvdr_sparse_basis` and `ssurvdr_basis`.
#'   It contains the row-sparse orthonormal basis, selected predictor indices,
#'   row-norm ranking, retained basis energy, and subspace diagnostics.
#'
#' @export
sparsify_sdr_basis <- function(
    basis_fit,
    support_size,
    variable_names = NULL,
    rank_tol = 1e-10
) {
  if (missing(basis_fit)) {
    stop("'basis_fit' must be supplied.", call. = FALSE)
  }

  if (!inherits(basis_fit, "ssurvdr_basis")) {
    stop(
      "'basis_fit' must inherit from 'ssurvdr_basis'.",
      call. = FALSE
    )
  }

  B <- as.matrix(basis_fit$basis)

  if (!is.numeric(B) || anyNA(B) || any(!is.finite(B))) {
    stop("The SDR basis must contain only finite numeric values.", call. = FALSE)
  }

  p <- nrow(B)
  d <- ncol(B)

  if (p < 1L || d < 1L) {
    stop("The SDR basis must have positive dimensions.", call. = FALSE)
  }

  if (!is.numeric(support_size) || length(support_size) != 1L) {
    stop("'support_size' must be a single positive integer.", call. = FALSE)
  }

  support_size <- as.integer(support_size)

  if (support_size < d) {
    stop(
      "'support_size' must be at least the structural dimension d.",
      call. = FALSE
    )
  }

  if (support_size > p) {
    stop("'support_size' cannot exceed the number of predictors.", call. = FALSE)
  }

  if (!is.numeric(rank_tol) || length(rank_tol) != 1L || rank_tol <= 0) {
    stop("'rank_tol' must be a single positive numeric value.", call. = FALSE)
  }

  if (is.null(variable_names)) {
    variable_names <- rownames(B)

    if (is.null(variable_names)) {
      variable_names <- paste0("V", seq_len(p))
    }
  }

  if (length(variable_names) != p) {
    stop(
      "'variable_names' must have length equal to the number of predictor rows.",
      call. = FALSE
    )
  }

  row_norms <- sqrt(rowSums(B^2))
  ranking <- order(row_norms, decreasing = TRUE)
  selected_indices <- ranking[seq_len(support_size)]

  total_energy <- sum(B^2)
  full_support <- support_size == p

  # The full-support endpoint is the identity element of the sparsification
  # path.  No rows are removed, so the original dense basis should be returned
  # exactly rather than being re-orthonormalised a second time.  This both
  # preserves the fitted dense model and avoids floating-point drift in the
  # benchmark endpoint used by the support-path acceptance tests.
  if (full_support) {
    B_sparse <- B

    retained_rank <- qr(
      B,
      tol = rank_tol
    )$rank

    if (retained_rank < d) {
      stop(
        paste0(
          "The full SDR basis has numerical rank ", retained_rank,
          ", which is smaller than d = ", d, "."
        ),
        call. = FALSE
      )
    }

    tolerance <- rank_tol
    orthogonality_error <- max(
      abs(
        crossprod(B_sparse) - diag(d)
      )
    )

    singular_values <- rep(1, d)
    projection_distance <- 0
    retained_energy <- total_energy
    retained_energy_proportion <- if (total_energy > 0) 1 else NA_real_
  } else {
    B_zero <- matrix(
      0,
      nrow = p,
      ncol = d,
      dimnames = dimnames(B)
    )

    B_zero[selected_indices, ] <- B[selected_indices, , drop = FALSE]

    gram <- crossprod(B_zero)
    gram <- (gram + t(gram)) / 2

    eig <- eigen(
      gram,
      symmetric = TRUE
    )

    maximum_eigenvalue <- max(eig$values, 0)
    tolerance <- max(
      rank_tol,
      maximum_eigenvalue * sqrt(.Machine$double.eps)
    )

    retained_rank <- sum(eig$values > tolerance)

    if (retained_rank < d) {
      stop(
        paste0(
          "The retained support of size ", support_size,
          " has numerical rank ", retained_rank,
          ", which is smaller than d = ", d, "."
        ),
        call. = FALSE
      )
    }

    inv_sqrt <- eig$vectors %*%
      diag(1 / sqrt(pmax(eig$values, tolerance)), nrow = d) %*%
      t(eig$vectors)

    B_sparse <- B_zero %*% inv_sqrt

    # Preserve exact row sparsity after floating-point matrix multiplication.
    excluded_indices <- setdiff(seq_len(p), selected_indices)

    if (length(excluded_indices) > 0L) {
      B_sparse[excluded_indices, ] <- 0
    }

    orthogonality_error <- max(
      abs(
        crossprod(B_sparse) - diag(d)
      )
    )

    # Orthonormalise the dense basis only for subspace diagnostics. This does not
    # alter the sparse basis used for modelling.
    dense_q <- qr.Q(
      qr(B),
      complete = FALSE
    )

    cross_basis <- crossprod(
      dense_q,
      B_sparse
    )

    singular_values <- svd(
      cross_basis,
      nu = 0,
      nv = 0
    )$d

    singular_values <- pmin(
      pmax(singular_values, 0),
      1
    )

    projection_distance <- sqrt(
      max(
        0,
        2 * d - 2 * sum(singular_values^2)
      )
    )

    retained_energy <- sum(B[selected_indices, , drop = FALSE]^2)

    retained_energy_proportion <- if (total_energy > 0) {
      retained_energy / total_energy
    } else {
      NA_real_
    }
  }

  score_table <- data.frame(
    variable = variable_names,
    index = seq_len(p),
    row_norm = row_norms,
    rank = rank(
      -row_norms,
      ties.method = "first"
    ),
    selected = seq_len(p) %in% selected_indices,
    stringsAsFactors = FALSE
  )

  score_table <- score_table[
    order(score_table$rank),
    ,
    drop = FALSE
  ]

  rownames(score_table) <- NULL

  selected_table <- score_table[
    score_table$selected,
    ,
    drop = FALSE
  ]

  out <- list(
    basis = B_sparse,
    projection_matrix = NULL,
    eigenvalues = basis_fit$eigenvalues,
    full_eigenvalues = basis_fit$full_eigenvalues,
    explained_ratio = basis_fit$explained_ratio,
    cumulative_ratio = basis_fit$cumulative_ratio,
    diagnostics = basis_fit$diagnostics,
    d = d,
    p = p,
    covariance_method = basis_fit$covariance_method,
    preprocessing = basis_fit$preprocessing,
    backend = basis_fit$backend,
    data_span_restricted = basis_fit$data_span_restricted,
    row_space_rank = basis_fit$row_space_rank,
    null_eigenvalue = basis_fit$null_eigenvalue,
    null_multiplicity = basis_fit$null_multiplicity,
    n_row_space_eigenvalues_above_null =
      basis_fit$n_row_space_eigenvalues_above_null,
    projection_stored = FALSE,
    sparse = TRUE,
    full_support_identity = full_support,
    support_size = support_size,
    selected_indices = selected_indices,
    selected_variables = variable_names[selected_indices],
    variable_scores = score_table,
    selected_scores = selected_table,
    retained_rank = retained_rank,
    rank_tolerance = tolerance,
    retained_basis_energy = retained_energy,
    total_basis_energy = total_energy,
    retained_energy_proportion = retained_energy_proportion,
    orthogonality_error = orthogonality_error,
    principal_cosines = singular_values,
    projection_distance_from_dense = projection_distance,
    source_support_size = p,
    call = match.call()
  )

  class(out) <- c(
    "ssurvdr_sparse_basis",
    "ssurvdr_basis"
  )

  out
}
