#' Sparse projection from SAVE kernel
#'
#' Computes sparse sufficient dimension reduction directions from a SAVE kernel
#' using weighted soft-thresholding followed by eigendecomposition.
#'
#' This function thresholds entries of a primal SAVE kernel before extracting
#' directions. It requires a dense primal kernel. For a row-support path using
#' a fitted primal or dual basis, use [fit_ssurvdr_sparse_path()].
#'
#' @param save_fit Object returned by [save_kernel()].
#' @param d Structural dimension.
#' @param lambda Non-negative sparsity parameter.
#' @param weights Optional non-negative variable weights of length p.
#' @param threshold Final row-norm threshold for exact zeroing of weak rows.
#' @param standardise Logical. If TRUE, normalise sparse directions to unit norm.
#'
#' @return A list containing sparse basis, sparse projection matrix, selected
#'   variables, thresholded kernel, and diagnostics.
#'
#' @export
sparse_projection <- function(
    save_fit,
    d = 1,
    lambda = 0.01,
    weights = NULL,
    threshold = 1e-6,
    standardise = TRUE
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

  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0) {
    stop("'lambda' must be a single non-negative numeric value.", call. = FALSE)
  }

  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0) {
    stop("'threshold' must be a single non-negative numeric value.", call. = FALSE)
  }

  d <- as.integer(d)

  if (identical(save_fit$backend, "dual") || is.null(save_fit$save_matrix)) {
    stop(
      paste0(
        "sparse_projection() currently requires a dense primal SAVE matrix. ",
        "A mathematically validated sparse operator for the compact dual ",
        "backend has not yet been implemented."
      ),
      call. = FALSE
    )
  }

  M <- save_fit$save_matrix

  if (!is.matrix(M) || !is.numeric(M)) {
    stop("SAVE matrix is invalid.", call. = FALSE)
  }

  if (nrow(M) != ncol(M)) {
    stop("SAVE matrix must be square.", call. = FALSE)
  }

  p <- ncol(M)

  if (d > p) {
    stop("'d' cannot exceed the number of predictors.", call. = FALSE)
  }

  if (is.null(weights)) {
    weights <- rep(1, p)
  }

  if (!is.numeric(weights) || length(weights) != p) {
    stop("'weights' must be NULL or a numeric vector of length p.", call. = FALSE)
  }

  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("'weights' must contain only finite non-missing values.", call. = FALSE)
  }

  if (any(weights < 0)) {
    stop("'weights' must contain non-negative values.", call. = FALSE)
  }

  M <- (M + t(M)) / 2

  penalty_matrix <- sqrt(outer(weights, weights))

  soft_threshold <- function(x, penalty) {
    sign(x) * pmax(abs(x) - penalty, 0)
  }

  M_sparse <- soft_threshold(
    x = M,
    penalty = lambda * penalty_matrix
  )

  M_sparse <- (M_sparse + t(M_sparse)) / 2

  diag(M_sparse) <- pmax(diag(M_sparse), 0)

  if (all(abs(M_sparse) <= .Machine$double.eps)) {
    warning(
      "All entries were thresholded to zero. Returning zero sparse basis.",
      call. = FALSE
    )

    basis <- matrix(0, nrow = p, ncol = d)
    projection <- matrix(0, nrow = p, ncol = p)
    eigenvalues <- rep(0, d)
  } else {
    eig <- eigen(
      M_sparse,
      symmetric = TRUE
    )

    ord <- order(eig$values, decreasing = TRUE)

    eigenvalues_all <- eig$values[ord]
    eigenvectors_all <- eig$vectors[, ord, drop = FALSE]

    basis <- eigenvectors_all[, seq_len(d), drop = FALSE]
    eigenvalues <- eigenvalues_all[seq_len(d)]

    row_norms_initial <- sqrt(rowSums(basis^2))

    weak_rows <- row_norms_initial <= threshold

    if (any(weak_rows)) {
      basis[weak_rows, ] <- 0
    }

    if (isTRUE(standardise)) {
      for (j in seq_len(ncol(basis))) {
        norm_j <- sqrt(sum(basis[, j]^2))

        if (norm_j > 0) {
          basis[, j] <- basis[, j] / norm_j
        }
      }
    }

    projection <- basis %*% t(basis)
  }

  row_norms <- sqrt(rowSums(basis^2))

  selected_indices <- which(row_norms > threshold)

  variable_scores <- data.frame(
    index = seq_len(p),
    score = row_norms,
    selected = row_norms > threshold,
    rank = rank(
      -row_norms,
      ties.method = "first"
    ),
    stringsAsFactors = FALSE
  )

  variable_scores <- variable_scores[
    order(variable_scores$rank),
    ,
    drop = FALSE
  ]

  rownames(variable_scores) <- NULL

  out <- list(
    basis = basis,
    projection_matrix = projection,
    sparse_kernel = M_sparse,
    original_kernel = M,
    eigenvalues = eigenvalues,
    variable_scores = variable_scores,
    selected_indices = selected_indices,
    n_selected = length(selected_indices),
    lambda = lambda,
    weights = weights,
    threshold = threshold,
    d = d,
    p = p,
    covariance_method = save_fit$covariance_method,
    save_fit = save_fit
  )

  class(out) <- "ssurvdr_sparse_projection"

  out
}
