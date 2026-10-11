#' Maximum Entropy Covariance Estimator
#'
#' Computes the Maximum Entropy Covariance Estimator (MEC) from a predictor
#' matrix and survival-informed slice membership. The estimator combines a
#' pooled within-slice covariance matrix with the entropy-maximising slice
#' covariance matrix, with optional diagonal stabilisation.
#'
#' @param X Numeric matrix or data frame of predictors with observations in rows.
#' @param slices Either an object returned by [survival_slicing()] or an integer
#'   vector of slice memberships.
#' @param rho Convex mixing weight for the entropy-selected covariance matrix.
#'   Must lie between 0 and 1, inclusive.
#' @param alpha Optional diagonal stabilisation weight. Must lie between 0 and 1, inclusive.
#' @param eps Small positive constant used for eigenvalue and logarithmic
#'   stabilisation.
#' @param centre Logical. If TRUE, centre columns of X before covariance
#'   estimation.
#' @param scale Logical. If TRUE, scale columns of X before covariance
#'   estimation.
#'
#' @return A list containing the MEC covariance matrix, pooled covariance matrix,
#'   entropy-selected covariance matrix, entropy scores, selected slice, and
#'   diagnostic quantities.
#'
#' @export
mec_covariance <- function(
    X,
    slices,
    rho = 0.5,
    alpha = 0.05,
    eps = 1e-8,
    centre = TRUE,
    scale = FALSE
) {
  if (missing(X)) {
    stop("'X' must be supplied.", call. = FALSE)
  }

  if (missing(slices)) {
    stop("'slices' must be supplied.", call. = FALSE)
  }

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    stop("'X' must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (!is.numeric(X)) {
    stop("'X' must contain numeric values only.", call. = FALSE)
  }

  if (anyNA(X)) {
    stop("'X' must not contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(X))) {
    stop("'X' must contain only finite values.", call. = FALSE)
  }

  if (!is.numeric(rho) || length(rho) != 1L || rho < 0 || rho > 1) {
    stop("'rho' must be a single numeric value in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("'eps' must be a single positive numeric value.", call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)

  if (n < 2L) {
    stop("'X' must contain at least two observations.", call. = FALSE)
  }

  if (p < 1L) {
    stop("'X' must contain at least one predictor.", call. = FALSE)
  }

  if (inherits(slices, "ssurvdr_slices")) {
    slice_id <- slices$slice_id
  } else {
    slice_id <- slices
  }

  if (length(slice_id) != n) {
    stop("'slices' must have length equal to nrow(X).", call. = FALSE)
  }

  if (anyNA(slice_id)) {
    stop("'slices' must not contain missing values.", call. = FALSE)
  }

  slice_id <- as.integer(as.factor(slice_id))
  unique_slices <- sort(unique(slice_id))

  if (length(unique_slices) < 1L) {
    stop("At least one slice is required.", call. = FALSE)
  }

  X_work <- X

  if (isTRUE(centre) || isTRUE(scale)) {
    X_work <- scale(
      X_work,
      center = isTRUE(centre),
      scale = isTRUE(scale)
    )
    X_work <- as.matrix(X_work)
  }

  slice_covariances <- vector("list", length(unique_slices))
  entropy_scores <- numeric(length(unique_slices))
  slice_sizes <- integer(length(unique_slices))
  slice_labels <- character(length(unique_slices))

  pooled_numerator <- matrix(0, nrow = p, ncol = p)
  pooled_denominator <- 0L

  stable_cov <- function(Z) {
    if (nrow(Z) <= 1L) {
      return(matrix(0, nrow = ncol(Z), ncol = ncol(Z)))
    }

    stats::cov(Z)
  }

  stable_entropy <- function(S) {
    S <- (S + t(S)) / 2

    eigen_values <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    eigen_values <- pmax(eigen_values, 0)

    sum(log(eigen_values + eps))
  }

  stable_condition_number <- function(S) {
    S <- (S + t(S)) / 2

    eigen_values <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    eigen_values <- pmax(eigen_values, 0)

    max_eig <- max(eigen_values)
    min_eig <- min(eigen_values)

    if (min_eig <= eps) {
      return(Inf)
    }

    max_eig / min_eig
  }

  for (i in seq_along(unique_slices)) {
    current_slice <- unique_slices[i]
    idx <- which(slice_id == current_slice)
    Z <- X_work[idx, , drop = FALSE]

    S_h <- stable_cov(Z)
    S_h <- (S_h + t(S_h)) / 2

    slice_covariances[[i]] <- S_h
    entropy_scores[i] <- stable_entropy(S_h)
    slice_sizes[i] <- length(idx)
    slice_labels[i] <- paste0("slice_", current_slice)

    if (length(idx) > 1L) {
      pooled_numerator <- pooled_numerator + (length(idx) - 1L) * S_h
      pooled_denominator <- pooled_denominator + (length(idx) - 1L)
    }
  }

  if (pooled_denominator > 0L) {
    Sigma_pooled <- pooled_numerator / pooled_denominator
  } else {
    Sigma_pooled <- stats::cov(X_work)
  }

  Sigma_pooled <- (Sigma_pooled + t(Sigma_pooled)) / 2

  selected_position <- which.max(entropy_scores)
  selected_slice <- unique_slices[selected_position]
  Sigma_mec_selected <- slice_covariances[[selected_position]]

  Sigma_mec <- rho * Sigma_mec_selected + (1 - rho) * Sigma_pooled
  Sigma_mec <- (Sigma_mec + t(Sigma_mec)) / 2

  diagonal_target <- (sum(diag(Sigma_mec)) / p) * diag(p)

  Sigma_mec_stabilised <- (1 - alpha) * Sigma_mec + alpha * diagonal_target
  Sigma_mec_stabilised <- (Sigma_mec_stabilised + t(Sigma_mec_stabilised)) / 2

  eigen_raw <- eigen(Sigma_mec, symmetric = TRUE, only.values = TRUE)$values
  eigen_stabilised <- eigen(
    Sigma_mec_stabilised,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  diagnostics <- data.frame(
    slice_id = unique_slices,
    slice_label = slice_labels,
    n = slice_sizes,
    entropy = entropy_scores,
    stringsAsFactors = FALSE
  )

  out <- list(
    covariance = Sigma_mec_stabilised,
    covariance_unstabilised = Sigma_mec,
    pooled_covariance = Sigma_pooled,
    entropy_selected_covariance = Sigma_mec_selected,
    entropy_scores = diagnostics,
    selected_slice_id = selected_slice,
    selected_slice_position = selected_position,
    rho = rho,
    alpha = alpha,
    eps = eps,
    eigenvalues_unstabilised = eigen_raw,
    eigenvalues = eigen_stabilised,
    condition_number_unstabilised = stable_condition_number(Sigma_mec),
    condition_number = stable_condition_number(Sigma_mec_stabilised),
    n = n,
    p = p,
    centre = centre,
    scale = scale
  )

  class(out) <- "ssurvdr_mec_covariance"

  out
}
