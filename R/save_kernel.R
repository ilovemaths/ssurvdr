#' Compute survival SAVE kernel matrix
#'
#' Computes the covariance-stabilised survival SAVE operator for sufficient
#' dimension reduction in censored survival data. For very high-dimensional
#' problems, a compact dual-space backend avoids materialising dense p x p
#' covariance and SAVE matrices and restricts estimation to the observed data
#' span.
#'
#' @param X Numeric predictor matrix with observations in rows.
#' @param slices Either an object returned by [survival_slicing()] or an
#'   integer vector of slice memberships.
#' @param covariance Character string specifying covariance estimator.
#'   Currently supports "sample" and "mec".
#' @param rho MEC shrinkage parameter.
#' @param alpha MEC diagonal stabilisation parameter.
#' @param eps Small positive stabilisation constant.
#' @param centre Logical. If TRUE, centre predictors.
#' @param scale Logical. If TRUE, scale predictors.
#' @param backend Computational backend. "primal" reproduces the dense p x p
#'   implementation. "dual" uses a compact row-space representation. "auto"
#'   uses the dual backend when p > max(n, 1000) and the primal backend
#'   otherwise.
#'
#' @return A list containing the SAVE representation, eigensystem, preprocessing
#'   metadata, and diagnostics.
#'
#' @export
save_kernel <- function(
    X,
    slices,
    covariance = c("sample", "mec"),
    rho = 0.5,
    alpha = 0.05,
    eps = 1e-8,
    centre = TRUE,
    scale = FALSE,
    backend = c("auto", "primal", "dual")
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
    stop("'X' must be a numeric matrix.", call. = FALSE)
  }

  if (!is.numeric(X)) {
    stop("'X' must contain numeric values only.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' must contain only finite non-missing values.", call. = FALSE)
  }

  covariance <- match.arg(covariance)
  backend <- match.arg(backend)

  n <- nrow(X)
  p <- ncol(X)

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

  if (backend == "auto") {
    backend <- if (p > max(n, 1000L)) "dual" else "primal"
  }

  prep <- .prepare_ssurvdr_predictors(
    X = X,
    centre = centre,
    scale = scale
  )

  X_work <- prep$X
  preprocessing <- prep$preprocessing

  unique_slices <- sort(unique(slice_id))
  H <- length(unique_slices)

  if (H < 1L) {
    stop("At least one slice is required.", call. = FALSE)
  }

  # -----------------------------------------------------------------------
  # Compact dual-space backend
  # -----------------------------------------------------------------------
  if (backend == "dual") {
    compact <- if (covariance == "sample") {
      .dual_sample_save(
        X_work = X_work,
        slice_id = slice_id,
        eps = eps
      )
    } else {
      .dual_mec_save(
        X_work = X_work,
        slice_id = slice_id,
        rho = rho,
        alpha = alpha,
        eps = eps
      )
    }

    out <- list(
      save_matrix = NULL,
      reduced_save_matrix = compact$reduced_save_matrix,
      eigenvalues = compact$eigenvalues,
      eigenvectors = compact$eigenvectors,
      row_space_basis = compact$row_space_basis,
      row_space_rank = compact$row_space_rank,
      null_eigenvalue = compact$null_eigenvalue,
      null_multiplicity = compact$null_multiplicity,
      data_span_restricted = compact$data_span_restricted,
      slice_covariances = compact$slice_covariances,
      slice_covariance_representation = compact$slice_covariance_representation,
      slice_weights = compact$slice_weights,
      covariance_method = covariance,
      diagnostics = compact$diagnostics,
      preprocessing = preprocessing,
      backend = backend,
      n = n,
      p = p,
      n_slices = H
    )

    if (covariance == "mec") {
      out$mec_fit <- compact$mec_fit
    }

    class(out) <- "ssurvdr_save_kernel"
    return(out)
  }

  # -----------------------------------------------------------------------
  # Dense primal backend
  # -----------------------------------------------------------------------
  slice_covariances <- vector("list", H)
  slice_weights <- numeric(H)

  stable_cov <- function(Z) {
    if (nrow(Z) <= 1L) {
      return(matrix(0, nrow = ncol(Z), ncol = ncol(Z)))
    }

    stats::cov(Z)
  }

  if (covariance == "mec") {
    mec_fit <- mec_covariance(
      X = X_work,
      slices = slice_id,
      rho = rho,
      alpha = alpha,
      eps = eps,
      centre = FALSE,
      scale = FALSE
    )

    global_covariance <- mec_fit$covariance
  }

  M_save <- matrix(0, nrow = p, ncol = p)

  diagnostics <- data.frame(
    slice_id = integer(H),
    n = integer(H),
    weight = numeric(H),
    condition_number = numeric(H),
    stringsAsFactors = FALSE
  )

  compute_condition_number <- function(S) {
    eigen_values <- eigen(
      (S + t(S)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    eigen_values <- pmax(eigen_values, eps)
    max(eigen_values) / min(eigen_values)
  }

  I_p <- diag(p)

  for (h in seq_along(unique_slices)) {
    current_slice <- unique_slices[h]
    idx <- which(slice_id == current_slice)
    Z_h <- X_work[idx, , drop = FALSE]
    p_h <- length(idx) / n

    if (covariance == "sample") {
      Sigma_h <- stable_cov(Z_h)
    } else {
      Sigma_h <- global_covariance
    }

    Sigma_h <- (Sigma_h + t(Sigma_h)) / 2

    slice_covariances[[h]] <- Sigma_h
    slice_weights[h] <- p_h

    delta <- I_p - Sigma_h
    M_save <- M_save + p_h * (delta %*% delta)

    diagnostics[h, ] <- list(
      slice_id = current_slice,
      n = length(idx),
      weight = p_h,
      condition_number = compute_condition_number(Sigma_h)
    )
  }

  M_save <- (M_save + t(M_save)) / 2

  eig <- eigen(M_save, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)

  out <- list(
    save_matrix = M_save,
    reduced_save_matrix = NULL,
    eigenvalues = eig$values[ord],
    eigenvectors = eig$vectors[, ord, drop = FALSE],
    row_space_basis = NULL,
    row_space_rank = min(n, p),
    null_eigenvalue = NA_real_,
    null_multiplicity = 0L,
    data_span_restricted = FALSE,
    slice_covariances = slice_covariances,
    slice_covariance_representation = "dense_primal",
    slice_weights = slice_weights,
    covariance_method = covariance,
    diagnostics = diagnostics,
    preprocessing = preprocessing,
    backend = backend,
    n = n,
    p = p,
    n_slices = H
  )

  if (covariance == "mec") {
    out$mec_fit <- mec_fit
  }

  class(out) <- "ssurvdr_save_kernel"
  out
}
