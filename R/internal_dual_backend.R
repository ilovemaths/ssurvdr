# Internal helpers for low-rank / dual-space ssurvdr computation.
# These functions are intentionally not exported.

.prepare_ssurvdr_predictors <- function(X, centre = TRUE, scale = FALSE) {
  X_work <- X

  if (isTRUE(centre) || isTRUE(scale)) {
    X_work <- base::scale(
      X_work,
      center = isTRUE(centre),
      scale = isTRUE(scale)
    )
    X_work <- as.matrix(X_work)
  }

  list(
    X = X_work,
    preprocessing = list(
      centre_requested = isTRUE(centre),
      scale_requested = isTRUE(scale),
      center = attr(X_work, "scaled:center"),
      scale = attr(X_work, "scaled:scale")
    )
  )
}

.apply_ssurvdr_preprocessing <- function(X, preprocessing = NULL) {
  X_work <- X

  if (is.null(preprocessing)) {
    return(X_work)
  }

  center_values <- preprocessing$center
  scale_values <- preprocessing$scale

  if (!is.null(center_values)) {
    X_work <- sweep(X_work, 2, center_values, FUN = "-")
  }

  if (!is.null(scale_values)) {
    X_work <- sweep(X_work, 2, scale_values, FUN = "/")
  }

  as.matrix(X_work)
}

.center_slice_matrix <- function(Z) {
  if (nrow(Z) <= 1L) {
    return(matrix(0, nrow = nrow(Z), ncol = ncol(Z)))
  }

  sweep(Z, 2, colMeans(Z), FUN = "-")
}

.row_space_basis <- function(A, eps = 1e-10) {
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("'A' must be a numeric matrix.", call. = FALSE)
  }

  p <- ncol(A)

  if (nrow(A) == 0L || p == 0L) {
    stop("'A' must have positive dimensions.", call. = FALSE)
  }

  gram <- tcrossprod(A)
  gram <- (gram + t(gram)) / 2

  eig <- eigen(gram, symmetric = TRUE)

  maximum <- if (length(eig$values) > 0L) {
    max(eig$values, 0)
  } else {
    0
  }

  tolerance <- max(eps, maximum * sqrt(.Machine$double.eps))
  keep <- which(eig$values > tolerance)

  if (length(keep) == 0L) {
    return(list(
      basis = matrix(0, nrow = p, ncol = 0L),
      eigenvalues = numeric(0),
      rank = 0L,
      tolerance = tolerance
    ))
  }

  U <- eig$vectors[, keep, drop = FALSE]
  values <- eig$values[keep]

  V <- crossprod(A, U)
  V <- sweep(V, 2, sqrt(values), FUN = "/")

  # A light numerical re-orthogonalisation in the small r-dimensional span.
  qr_v <- qr(V)
  V <- qr.Q(qr_v, complete = FALSE)

  list(
    basis = V,
    eigenvalues = values,
    rank = ncol(V),
    tolerance = tolerance
  )
}

.nonzero_cov_eigenvalues <- function(Y, denominator, eps = 1e-10) {
  if (nrow(Y) <= 1L || denominator <= 0) {
    return(numeric(0))
  }

  gram <- tcrossprod(Y) / denominator
  gram <- (gram + t(gram)) / 2

  values <- eigen(
    gram,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  maximum <- max(values, 0)
  tolerance <- max(eps, maximum * sqrt(.Machine$double.eps))

  sort(values[values > tolerance], decreasing = TRUE)
}

.compact_entropy <- function(nonzero_values, p, eps = 1e-8) {
  r <- length(nonzero_values)
  zero_count <- max(0L, p - r)

  sum(log(pmax(nonzero_values, 0) + eps)) +
    zero_count * log(eps)
}

.compact_condition_number <- function(nonzero_values, p, eps = 1e-8) {
  if (p < 1L) {
    return(NA_real_)
  }

  if (length(nonzero_values) == 0L) {
    return(1)
  }

  maximum <- max(nonzero_values, eps)

  # In p > rank settings, the dense implementation clamps zero eigenvalues
  # to eps. Preserve that diagnostic convention without materialising p x p.
  minimum <- if (length(nonzero_values) < p) {
    eps
  } else {
    max(min(nonzero_values), eps)
  }

  maximum / minimum
}

.dual_sample_save <- function(X_work, slice_id, eps = 1e-8) {
  n <- nrow(X_work)
  p <- ncol(X_work)
  unique_slices <- sort(unique(slice_id))
  H <- length(unique_slices)

  centered_slices <- vector("list", H)
  slice_weights <- numeric(H)
  diagnostics <- data.frame(
    slice_id = integer(H),
    n = integer(H),
    weight = numeric(H),
    condition_number = numeric(H),
    stringsAsFactors = FALSE
  )

  stacked_rows <- vector("list", H)

  for (h in seq_along(unique_slices)) {
    current_slice <- unique_slices[h]
    idx <- which(slice_id == current_slice)
    Z_h <- X_work[idx, , drop = FALSE]
    Y_h <- .center_slice_matrix(Z_h)

    centered_slices[[h]] <- Y_h
    stacked_rows[[h]] <- Y_h
    slice_weights[h] <- length(idx) / n

    nonzero_values <- .nonzero_cov_eigenvalues(
      Y_h,
      denominator = max(length(idx) - 1L, 1L),
      eps = eps
    )

    diagnostics[h, ] <- list(
      slice_id = current_slice,
      n = length(idx),
      weight = slice_weights[h],
      condition_number = .compact_condition_number(
        nonzero_values,
        p = p,
        eps = eps
      )
    )
  }

  A <- do.call(rbind, stacked_rows)
  span_fit <- .row_space_basis(A, eps = eps)
  V <- span_fit$basis
  r <- span_fit$rank

  if (r == 0L) {
    stop("The slice-centred predictor matrix has zero numerical rank.", call. = FALSE)
  }

  M_reduced <- diag(r)
  M_reduced[,] <- 0

  for (h in seq_along(unique_slices)) {
    Y_h <- centered_slices[[h]]
    m_h <- nrow(Y_h)

    if (m_h <= 1L) {
      S_h_reduced <- matrix(0, nrow = r, ncol = r)
    } else {
      YV <- Y_h %*% V
      S_h_reduced <- crossprod(YV) / (m_h - 1L)
      S_h_reduced <- (S_h_reduced + t(S_h_reduced)) / 2
    }

    delta <- diag(r) - S_h_reduced
    M_reduced <- M_reduced +
      slice_weights[h] * (delta %*% delta)
  }

  M_reduced <- (M_reduced + t(M_reduced)) / 2

  eig <- eigen(M_reduced, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  values <- eig$values[ord]
  vectors_reduced <- eig$vectors[, ord, drop = FALSE]
  vectors_predictor <- V %*% vectors_reduced

  list(
    reduced_save_matrix = M_reduced,
    eigenvalues = values,
    eigenvectors = vectors_predictor,
    row_space_basis = V,
    row_space_rank = r,
    null_eigenvalue = 1,
    null_multiplicity = max(0L, p - r),
    slice_covariances = NULL,
    slice_covariance_representation = "compact_row_space",
    slice_weights = slice_weights,
    diagnostics = diagnostics,
    data_span_restricted = p > r
  )
}

.dual_mec_save <- function(
    X_work,
    slice_id,
    rho = 0.5,
    alpha = 0.05,
    eps = 1e-8
) {
  n <- nrow(X_work)
  p <- ncol(X_work)
  unique_slices <- sort(unique(slice_id))
  H <- length(unique_slices)

  centered_slices <- vector("list", H)
  entropy_scores <- numeric(H)
  slice_sizes <- integer(H)
  slice_weights <- numeric(H)
  diagnostics <- data.frame(
    slice_id = integer(H),
    n = integer(H),
    weight = numeric(H),
    condition_number = numeric(H),
    stringsAsFactors = FALSE
  )

  pooled_rows <- list()
  pooled_position <- 0L
  pooled_denominator <- 0L

  for (h in seq_along(unique_slices)) {
    current_slice <- unique_slices[h]
    idx <- which(slice_id == current_slice)
    Z_h <- X_work[idx, , drop = FALSE]
    Y_h <- .center_slice_matrix(Z_h)
    m_h <- length(idx)

    centered_slices[[h]] <- Y_h
    slice_sizes[h] <- m_h
    slice_weights[h] <- m_h / n

    nonzero_values <- .nonzero_cov_eigenvalues(
      Y_h,
      denominator = max(m_h - 1L, 1L),
      eps = eps
    )

    entropy_scores[h] <- .compact_entropy(
      nonzero_values,
      p = p,
      eps = eps
    )

    diagnostics[h, ] <- list(
      slice_id = current_slice,
      n = m_h,
      weight = slice_weights[h],
      condition_number = .compact_condition_number(
        nonzero_values,
        p = p,
        eps = eps
      )
    )

    if (m_h > 1L) {
      pooled_position <- pooled_position + 1L
      pooled_rows[[pooled_position]] <- Y_h
      pooled_denominator <- pooled_denominator + (m_h - 1L)
    }
  }

  if (pooled_denominator <= 0L) {
    stop("Unable to construct pooled within-slice covariance.", call. = FALSE)
  }

  selected_position <- which.max(entropy_scores)
  selected_slice <- unique_slices[selected_position]
  Y_selected <- centered_slices[[selected_position]]
  m_selected <- nrow(Y_selected)

  pooled_stack <- do.call(rbind, pooled_rows)
  A_pool <- pooled_stack / sqrt(pooled_denominator)

  if (m_selected > 1L) {
    A_selected <- Y_selected / sqrt(m_selected - 1L)
  } else {
    A_selected <- matrix(0, nrow = 0L, ncol = p)
  }

  A_mix <- rbind(
    sqrt(rho) * A_selected,
    sqrt(1 - rho) * A_pool
  )

  span_fit <- .row_space_basis(A_mix, eps = eps)
  V <- span_fit$basis
  r <- span_fit$rank

  if (r == 0L) {
    stop("The MEC low-rank covariance component has zero numerical rank.", call. = FALSE)
  }

  AV <- A_mix %*% V
  S_reduced <- crossprod(AV)
  S_reduced <- (S_reduced + t(S_reduced)) / 2

  trace_low_rank <- sum(A_mix^2)
  tau <- trace_low_rank / p
  diagonal_shift <- alpha * tau
  low_rank_multiplier <- 1 - alpha

  Sigma_reduced <-
    low_rank_multiplier * S_reduced +
    diagonal_shift * diag(r)

  base_null_component <- 1 - diagonal_shift
  delta_reduced <-
    base_null_component * diag(r) -
    low_rank_multiplier * S_reduced

  M_reduced <- delta_reduced %*% delta_reduced
  M_reduced <- (M_reduced + t(M_reduced)) / 2

  eig <- eigen(M_reduced, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  values <- eig$values[ord]
  vectors_reduced <- eig$vectors[, ord, drop = FALSE]
  vectors_predictor <- V %*% vectors_reduced

  sigma_nonzero <- eigen(
    Sigma_reduced,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  null_sigma_value <- diagonal_shift
  max_sigma <- max(c(sigma_nonzero, null_sigma_value), na.rm = TRUE)
  min_sigma <- min(c(sigma_nonzero, null_sigma_value), na.rm = TRUE)

  condition_number <- if (min_sigma <= eps) {
    Inf
  } else {
    max_sigma / min_sigma
  }

  entropy_diagnostics <- data.frame(
    slice_id = unique_slices,
    slice_label = paste0("slice_", unique_slices),
    n = slice_sizes,
    entropy = entropy_scores,
    stringsAsFactors = FALSE
  )

  mec_fit <- list(
    covariance = NULL,
    covariance_unstabilised = NULL,
    pooled_covariance = NULL,
    entropy_selected_covariance = NULL,
    entropy_scores = entropy_diagnostics,
    selected_slice_id = selected_slice,
    selected_slice_position = selected_position,
    rho = rho,
    alpha = alpha,
    eps = eps,
    trace_target = tau,
    diagonal_shift = diagonal_shift,
    row_space_rank = r,
    condition_number = condition_number,
    representation = "compact_low_rank_plus_diagonal",
    n = n,
    p = p
  )
  class(mec_fit) <- "ssurvdr_mec_covariance"

  list(
    reduced_save_matrix = M_reduced,
    eigenvalues = values,
    eigenvectors = vectors_predictor,
    row_space_basis = V,
    row_space_rank = r,
    null_eigenvalue = base_null_component^2,
    null_multiplicity = max(0L, p - r),
    slice_covariances = NULL,
    slice_covariance_representation = "compact_low_rank_plus_diagonal",
    slice_weights = slice_weights,
    diagnostics = diagnostics,
    data_span_restricted = p > r,
    mec_fit = mec_fit
  )
}
