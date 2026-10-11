#' Extract sufficient dimension reduction basis
#'
#' Extracts leading sufficient dimension reduction directions from an object
#' returned by [save_kernel()]. The function reuses the eigensystem already
#' computed by save_kernel() and therefore avoids a redundant eigendecomposition.
#'
#' For a dual-space SAVE fit, estimation is restricted to the observed data span.
#' This removes non-identifiable null-space directions that arise when p > n.
#'
#' @param save_fit Object returned by [save_kernel()].
#' @param d Structural dimension to extract.
#' @param standardise Logical. If TRUE, normalise basis vectors to unit norm.
#' @param sort_eigenvalues Logical. If TRUE, sort eigenvalues in decreasing order.
#' @param store_projection Logical or NULL. If TRUE, store B B'. If NULL, dual
#'   fits never store the dense p x p projection matrix; primal fits store it only
#'   when p <= 2000. This prevents dual-space fits from reintroducing a dense
#'   high-dimensional allocation after basis extraction.
#'
#' @return A list containing SDR basis vectors, eigenvalues, optional projection
#'   matrix, preprocessing metadata, and diagnostics.
#'
#' @export
extract_sdr_basis <- function(
    save_fit,
    d = 1,
    standardise = TRUE,
    sort_eigenvalues = TRUE,
    store_projection = NULL
) {
  if (missing(save_fit)) {
    stop("'save_fit' must be supplied.", call. = FALSE)
  }

  if (!inherits(save_fit, "ssurvdr_save_kernel")) {
    stop(
      "'save_fit' must be an object returned by save_kernel().",
      call. = FALSE
    )
  }

  if (!is.numeric(d) || length(d) != 1L || d < 1) {
    stop("'d' must be a positive integer.", call. = FALSE)
  }

  d <- as.integer(d)
  p <- save_fit$p

  if (is.null(p)) {
    if (!is.null(save_fit$save_matrix)) {
      p <- ncol(save_fit$save_matrix)
    } else if (!is.null(save_fit$eigenvectors)) {
      p <- nrow(save_fit$eigenvectors)
    } else {
      stop("Could not determine the number of predictors.", call. = FALSE)
    }
  }

  if (d > p) {
    stop("'d' cannot exceed the number of predictors.", call. = FALSE)
  }

  eigenvalues <- save_fit$eigenvalues
  eigenvectors <- save_fit$eigenvectors

  # Backward compatibility for legacy save_kernel objects.
  if (is.null(eigenvalues) || is.null(eigenvectors)) {
    M <- save_fit$save_matrix

    if (!is.matrix(M)) {
      stop("SAVE eigensystem is unavailable and SAVE matrix is invalid.", call. = FALSE)
    }

    eig <- eigen(M, symmetric = TRUE)
    eigenvalues <- eig$values
    eigenvectors <- eig$vectors
  }

  if (!is.numeric(eigenvalues) || !is.matrix(eigenvectors)) {
    stop("SAVE eigensystem is invalid.", call. = FALSE)
  }

  if (nrow(eigenvectors) != p) {
    stop("SAVE eigenvectors have incompatible predictor dimension.", call. = FALSE)
  }

  if (isTRUE(sort_eigenvalues)) {
    ord <- order(eigenvalues, decreasing = TRUE)
    eigenvalues <- eigenvalues[ord]
    eigenvectors <- eigenvectors[, ord, drop = FALSE]
  }

  available_directions <- ncol(eigenvectors)

  if (d > available_directions) {
    stop(
      paste0(
        "Requested d = ", d,
        " exceeds the ", available_directions,
        " identifiable directions available from the ",
        if (identical(save_fit$backend, "dual")) "data-span-restricted dual" else "SAVE",
        " eigensystem."
      ),
      call. = FALSE
    )
  }

  selected_values <- eigenvalues[seq_len(d)]
  basis <- eigenvectors[, seq_len(d), drop = FALSE]

  if (isTRUE(standardise)) {
    for (j in seq_len(ncol(basis))) {
      norm_j <- sqrt(sum(basis[, j]^2))

      if (norm_j > 0) {
        basis[, j] <- basis[, j] / norm_j
      }
    }
  }

  if (is.null(store_projection)) {
    store_projection <- (
      !identical(save_fit$backend, "dual") &&
        p <= 2000L
    )
  }

  if (!is.logical(store_projection) || length(store_projection) != 1L) {
    stop("'store_projection' must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  projection_matrix <- if (isTRUE(store_projection)) {
    basis %*% t(basis)
  } else {
    NULL
  }

  positive_sum <- sum(pmax(eigenvalues, 0))

  if (positive_sum > 0) {
    explained_all <- pmax(eigenvalues, 0) / positive_sum
  } else {
    explained_all <- rep(0, length(eigenvalues))
  }

  explained_ratio <- explained_all[seq_len(d)]
  cumulative_ratio <- cumsum(explained_ratio)

  diagnostics <- data.frame(
    component = seq_len(length(eigenvalues)),
    eigenvalue = eigenvalues,
    explained_ratio = explained_all,
    cumulative_ratio = cumsum(explained_all),
    stringsAsFactors = FALSE
  )

  null_eigenvalue <- save_fit$null_eigenvalue
  null_multiplicity <- save_fit$null_multiplicity
  data_span_restricted <- isTRUE(save_fit$data_span_restricted)

  n_row_space_above_null <- NA_integer_

  if (
    data_span_restricted &&
    is.numeric(null_eigenvalue) &&
    length(null_eigenvalue) == 1L &&
    is.finite(null_eigenvalue)
  ) {
    tolerance <- max(1e-10, abs(null_eigenvalue) * sqrt(.Machine$double.eps))
    n_row_space_above_null <- sum(eigenvalues > null_eigenvalue + tolerance)

    if (d > n_row_space_above_null) {
      warning(
        paste0(
          "The requested structural dimension includes data-span directions ",
          "whose SAVE eigenvalues do not exceed the null-space eigenvalue. ",
          "The dual backend deliberately excludes arbitrary null-space ",
          "directions because they contain no observed predictor variation."
        ),
        call. = FALSE
      )
    }
  }

  out <- list(
    basis = basis,
    projection_matrix = projection_matrix,
    eigenvalues = selected_values,
    full_eigenvalues = eigenvalues,
    explained_ratio = explained_ratio,
    cumulative_ratio = cumulative_ratio,
    diagnostics = diagnostics,
    d = d,
    p = p,
    covariance_method = save_fit$covariance_method,
    preprocessing = save_fit$preprocessing,
    backend = save_fit$backend,
    data_span_restricted = data_span_restricted,
    row_space_rank = save_fit$row_space_rank,
    null_eigenvalue = null_eigenvalue,
    null_multiplicity = null_multiplicity,
    n_row_space_eigenvalues_above_null = n_row_space_above_null,
    projection_stored = isTRUE(store_projection)
  )

  class(out) <- "ssurvdr_basis"
  out
}
