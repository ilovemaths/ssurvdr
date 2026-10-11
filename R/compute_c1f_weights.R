#' Compute basis-adaptive variable weights
#'
#' Computes adaptive variable weights using row-wise initial SDR signal strength
#' and covariance-complexity adjustment. The function is intended for weighted
#' sparse projection in high-dimensional sufficient dimension reduction.
#'
#' @param basis_fit Object returned by [extract_sdr_basis()], an
#'   ssurvdr_sparse_projection object, or a numeric basis matrix.
#' @param matrix Optional covariance, SAVE, or information matrix used to compute
#'   eigenvalue complexity.
#' @param alpha Non-negative complexity adjustment parameter.
#' @param gamma Positive exponent applied to the inverse initial signal.
#' @param eps Small positive stabilisation constant.
#' @param normalise Logical. If TRUE, weights are normalised to have mean 1.
#' @param cap Optional positive upper cap for weights.
#'
#' @return A list containing adaptive weights, initial signal scores,
#'   complexity scores, and diagnostics.
#'
#' @details The returned `c1f` field is retained for compatibility. Its implemented
#' formula is `(p / 2) * log(arithmetic_mean / geometric_mean)` of the stabilised
#' eigenvalues, the same C1 expression used by [compute_icomp()]. The variable
#' weights use row norms and a separate loading-weighted eigenvalue-dispersion
#' score; the scalar `c1f` is a diagnostic and does not directly enter the weights.
#' With `matrix = NULL`, the complexity scores are zero and `c1f` is missing.
#' Finite weights are strictly positive, with no universal upper bound when
#' `cap = NULL`. Normalisation precedes optional capping.
#'
#' @examples
#' basis <- matrix(c(1, 2, 4), ncol = 1)
#' weights <- compute_c1f_weights(basis, matrix = diag(c(1, 2, 4)))
#' weights$weights
#' weights$c1f
#'
#' @export
compute_c1f_weights <- function(
    basis_fit,
    matrix = NULL,
    alpha = 1,
    gamma = 1,
    eps = 1e-8,
    normalise = TRUE,
    cap = NULL
) {
  if (missing(basis_fit)) {
    stop("'basis_fit' must be supplied.", call. = FALSE)
  }

  if (inherits(basis_fit, "ssurvdr_basis")) {
    B <- basis_fit$basis
  } else if (inherits(basis_fit, "ssurvdr_sparse_projection")) {
    B <- basis_fit$basis
  } else {
    B <- basis_fit
  }

  if (!is.matrix(B)) {
    stop("'basis_fit' must be an SDR basis object or a numeric matrix.", call. = FALSE)
  }

  if (!is.numeric(B)) {
    stop("The basis matrix must be numeric.", call. = FALSE)
  }

  if (anyNA(B) || any(!is.finite(B))) {
    stop("The basis matrix must contain only finite non-missing values.", call. = FALSE)
  }

  p <- nrow(B)

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0) {
    stop("'alpha' must be a single non-negative numeric value.", call. = FALSE)
  }

  if (!is.numeric(gamma) || length(gamma) != 1L || gamma <= 0) {
    stop("'gamma' must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("'eps' must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.null(cap)) {
    if (!is.numeric(cap) || length(cap) != 1L || cap <= 0) {
      stop("'cap' must be NULL or a single positive numeric value.", call. = FALSE)
    }
  }

  initial_signal <- sqrt(rowSums(B^2))

  if (is.null(matrix)) {
    complexity_score <- rep(0, p)

    complexity_source <- "none"
  } else {
    M <- as.matrix(matrix)

    if (!is.matrix(M) || !is.numeric(M)) {
      stop("'matrix' must be NULL or a numeric square matrix.", call. = FALSE)
    }

    if (nrow(M) != ncol(M)) {
      stop("'matrix' must be square.", call. = FALSE)
    }

    if (nrow(M) != p) {
      stop(
        "'matrix' dimension must match the number of rows in the basis matrix.",
        call. = FALSE
      )
    }

    if (anyNA(M) || any(!is.finite(M))) {
      stop("'matrix' must contain only finite non-missing values.", call. = FALSE)
    }

    M <- (M + t(M)) / 2

    eig <- eigen(M, symmetric = TRUE)

    eigenvalues <- pmax(eig$values, eps)
    eigenvectors <- eig$vectors

    mean_eigenvalue <- mean(eigenvalues)

    eigen_complexity <- (
      (eigenvalues - mean_eigenvalue)^2
    ) / (mean_eigenvalue^2 + eps)

    loading_energy <- eigenvectors^2

    complexity_score <- as.numeric(
      loading_energy %*% eigen_complexity
    )

    complexity_score <- pmax(complexity_score, 0)

    if (mean(complexity_score) > 0) {
      complexity_score <- complexity_score / mean(complexity_score)
    }

    complexity_source <- "matrix"
  }

  weights <- (1 + alpha * complexity_score) /
    ((initial_signal + eps)^gamma)

  if (isTRUE(normalise)) {
    weights <- weights / mean(weights)
  }

  if (!is.null(cap)) {
    weights <- pmin(weights, cap)
  }

  weights <- pmax(weights, eps)

  diagnostics <- data.frame(
    index = seq_len(p),
    initial_signal = initial_signal,
    complexity_score = complexity_score,
    weight = weights,
    rank_signal = rank(
      -initial_signal,
      ties.method = "first"
    ),
    rank_weight = rank(
      weights,
      ties.method = "first"
    ),
    stringsAsFactors = FALSE
  )

  diagnostics <- diagnostics[
    order(diagnostics$rank_weight),
    ,
    drop = FALSE
  ]

  rownames(diagnostics) <- NULL

  c1f_value <- if (!is.null(matrix)) {
    M <- (as.matrix(matrix) + t(as.matrix(matrix))) / 2
    eigenvalues <- pmax(
      eigen(M, symmetric = TRUE, only.values = TRUE)$values,
      eps
    )

    s <- length(eigenvalues)
    am <- mean(eigenvalues)
    gm <- exp(mean(log(eigenvalues)))

    (s / 2) * log(am / gm)
  } else {
    NA_real_
  }

  out <- list(
    weights = weights,
    initial_signal = initial_signal,
    complexity_score = complexity_score,
    diagnostics = diagnostics,
    c1f = c1f_value,
    alpha = alpha,
    gamma = gamma,
    eps = eps,
    normalise = normalise,
    cap = cap,
    complexity_source = complexity_source,
    p = p
  )

  class(out) <- "ssurvdr_c1f_weights"

  out
}
