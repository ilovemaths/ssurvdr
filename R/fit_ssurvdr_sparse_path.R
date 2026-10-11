#' Fit an ICOMP-guided row-sparsity path
#'
#' Evaluates a sequence of row-sparse SDR bases obtained from a previously fitted
#' dense SDR basis. Predictor rows are ranked by their Euclidean row norms, the
#' leading `k` rows are retained for each candidate support size, and the retained
#' basis is re-orthonormalised before refitting the reduced Cox model.
#'
#' Model comparison uses the existing Cox-model criteria returned by
#' [compute_icomp()]. In particular, `ICOMP_IFIM` and `ICOMP_MISSPEC` remain the
#' package's established information-complexity criteria. No additional support-
#' size penalty is silently added. Consequently, the returned path should be
#' inspected for boundary selection and compared with the dense-support candidate.
#'
#' This procedure is a post-estimation row-sparsification path. It is not a
#' group-lasso, Fantope, or penalised SAVE estimator.
#'
#' @param X Numeric predictor matrix with observations in rows.
#' @param time Numeric vector of observed survival times.
#' @param status Numeric or integer censoring indicator vector.
#' @param dense_fit Either an `ssurvdr_fit` object or an `ssurvdr_basis` object
#'   representing the dense SDR solution whose row norms define the support path.
#' @param support_grid Integer vector of candidate support sizes. If NULL, a
#'   deterministic default grid spanning from d to p is used.
#' @param selector Criterion used to rank the support path. One of
#'   `ICOMP_IFIM`, `ICOMP_MISSPEC`, `AIC`, `BIC`, or `CAIC`.
#' @param ties Character string specifying ties handling for the Cox model.
#' @param eps Small positive stabilisation constant passed to [compute_icomp()].
#' @param variable_names Optional character vector of predictor names.
#' @param rank_tol Positive numerical rank tolerance passed to
#'   `sparsify_sdr_basis()`.
#' @param continue_on_error Logical. If TRUE, invalid support sizes are recorded
#'   and the remaining path continues.
#' @param store_fits Logical. If TRUE, retain all successful sparse Cox fits.
#'   Otherwise only the selected fit and path table are returned.
#'
#' @return An object of class `ssurvdr_sparse_path` containing the support-path
#'   results, selected support, best sparse support excluding the fully dense
#'   candidate, and the selected fitted model.
#'
#' @export
fit_ssurvdr_sparse_path <- function(
    X,
    time,
    status,
    dense_fit,
    support_grid = NULL,
    selector = c(
      "ICOMP_IFIM",
      "ICOMP_MISSPEC",
      "AIC",
      "BIC",
      "CAIC"
    ),
    ties = "efron",
    eps = 1e-8,
    variable_names = NULL,
    rank_tol = 1e-10,
    continue_on_error = TRUE,
    store_fits = FALSE
) {
  if (missing(X) || missing(time) || missing(status) || missing(dense_fit)) {
    stop(
      "'X', 'time', 'status', and 'dense_fit' must all be supplied.",
      call. = FALSE
    )
  }

  selector <- match.arg(selector)

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("'X' must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' must contain only finite non-missing values.", call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)

  if (length(time) != n || length(status) != n) {
    stop("'time' and 'status' must have length nrow(X).", call. = FALSE)
  }

  if (inherits(dense_fit, "ssurvdr_fit")) {
    basis_fit <- dense_fit$basis_fit

    dense_reference <- list(
      criteria = dense_fit$criteria,
      selected = dense_fit$selected,
      diagnostics = dense_fit$diagnostics
    )
  } else if (inherits(dense_fit, "ssurvdr_basis")) {
    basis_fit <- dense_fit
    dense_reference <- NULL
  } else {
    stop(
      "'dense_fit' must inherit from 'ssurvdr_fit' or 'ssurvdr_basis'.",
      call. = FALSE
    )
  }

  if (!inherits(basis_fit, "ssurvdr_basis")) {
    stop("The dense fit does not contain a valid SDR basis.", call. = FALSE)
  }

  if (nrow(basis_fit$basis) != p) {
    stop(
      "The dense basis predictor dimension does not match ncol(X).",
      call. = FALSE
    )
  }

  d <- ncol(basis_fit$basis)

  if (is.null(variable_names)) {
    variable_names <- colnames(X)

    if (is.null(variable_names)) {
      variable_names <- rownames(basis_fit$basis)
    }

    if (is.null(variable_names)) {
      variable_names <- paste0("V", seq_len(p))
    }
  }

  if (length(variable_names) != p) {
    stop("'variable_names' must have length p.", call. = FALSE)
  }

  if (is.null(support_grid)) {
    support_grid <- unique(
      c(
        d,
        10L,
        20L,
        30L,
        40L,
        50L,
        75L,
        100L,
        150L,
        200L,
        300L,
        500L,
        750L,
        1000L,
        1500L,
        2000L,
        3000L,
        5000L,
        p
      )
    )
  }

  if (!is.numeric(support_grid) || length(support_grid) < 1L) {
    stop("'support_grid' must be a non-empty numeric/integer vector.", call. = FALSE)
  }

  support_grid <- sort(
    unique(
      as.integer(support_grid)
    )
  )

  support_grid <- support_grid[
    support_grid >= d & support_grid <= p
  ]

  if (length(support_grid) < 1L) {
    stop(
      "No support-grid values remain after enforcing d <= k <= p.",
      call. = FALSE
    )
  }

  # The dense candidate is required so that sparsification is compared against
  # the original full-support solution rather than against an arbitrary truncated
  # maximum support.
  if (!p %in% support_grid) {
    support_grid <- sort(unique(c(support_grid, p)))
  }

  if (!is.logical(continue_on_error) || length(continue_on_error) != 1L) {
    stop("'continue_on_error' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(store_fits) || length(store_fits) != 1L) {
    stop("'store_fits' must be TRUE or FALSE.", call. = FALSE)
  }

  results <- data.frame(
    support_size = support_grid,
    support_fraction = support_grid / p,
    retained_energy_proportion = NA_real_,
    projection_distance_from_dense = NA_real_,
    orthogonality_error = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    CAIC = NA_real_,
    ICOMP_IFIM = NA_real_,
    ICOMP_MISSPEC = NA_real_,
    c_index = NA_real_,
    c_index_se = NA_real_,
    runtime_seconds = NA_real_,
    is_dense_support = support_grid == p,
    status = "not_run",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  fit_cache <- vector("list", length(support_grid))

  for (i in seq_along(support_grid)) {
    k <- support_grid[i]
    start_i <- Sys.time()

    fit_i <- tryCatch(
      {
        sparse_basis <- sparsify_sdr_basis(
          basis_fit = basis_fit,
          support_size = k,
          variable_names = variable_names,
          rank_tol = rank_tol
        )

        cox_fit <- fit_reduced_cox(
          X = X,
          time = time,
          status = status,
          basis_fit = sparse_basis,
          ties = ties
        )

        icomp_fit <- compute_icomp(
          cox_fit = cox_fit,
          eps = eps
        )

        list(
          basis_fit = sparse_basis,
          cox_fit = cox_fit,
          icomp_fit = icomp_fit
        )
      },
      error = function(e) {
        e
      }
    )

    runtime_i <- as.numeric(
      difftime(
        Sys.time(),
        start_i,
        units = "secs"
      )
    )

    results$runtime_seconds[i] <- runtime_i

    if (inherits(fit_i, "error")) {
      results$status[i] <- "failed"
      results$error_message[i] <- conditionMessage(fit_i)

      if (!isTRUE(continue_on_error)) {
        stop(
          paste0(
            "Sparse support k = ", k,
            " failed: ", conditionMessage(fit_i)
          ),
          call. = FALSE
        )
      }

      next
    }

    criteria <- fit_i$icomp_fit$criteria

    criterion_value <- function(name) {
      value <- criteria$value[criteria$criterion == name]

      if (length(value) == 0L) {
        return(NA_real_)
      }

      as.numeric(value[1L])
    }

    concordance <- as.numeric(fit_i$cox_fit$concordance)

    results$retained_energy_proportion[i] <-
      fit_i$basis_fit$retained_energy_proportion
    results$projection_distance_from_dense[i] <-
      fit_i$basis_fit$projection_distance_from_dense
    results$orthogonality_error[i] <-
      fit_i$basis_fit$orthogonality_error
    results$AIC[i] <- criterion_value("AIC")
    results$BIC[i] <- criterion_value("BIC")
    results$CAIC[i] <- criterion_value("CAIC")
    results$ICOMP_IFIM[i] <- criterion_value("ICOMP_IFIM")
    results$ICOMP_MISSPEC[i] <- criterion_value("ICOMP_MISSPEC")
    results$c_index[i] <- if (length(concordance) >= 1L) concordance[1L] else NA_real_
    results$c_index_se[i] <- if (length(concordance) >= 2L) concordance[2L] else NA_real_
    results$status[i] <- "ok"
    results$error_message[i] <- NA_character_

    fit_cache[[i]] <- fit_i
  }

  successful <- which(
    results$status == "ok" &
      is.finite(results[[selector]])
  )

  if (length(successful) == 0L) {
    stop("No successful sparse-support model has a finite selector value.", call. = FALSE)
  }

  best_index <- successful[
    which.min(
      results[[selector]][successful]
    )
  ]

  best_row <- results[
    best_index,
    ,
    drop = FALSE
  ]

  sparse_successful <- successful[
    !results$is_dense_support[successful]
  ]

  if (length(sparse_successful) > 0L) {
    best_sparse_index <- sparse_successful[
      which.min(
        results[[selector]][sparse_successful]
      )
    ]

    best_sparse_row <- results[
      best_sparse_index,
      ,
      drop = FALSE
    ]
  } else {
    best_sparse_index <- NA_integer_
    best_sparse_row <- NULL
  }

  dense_index <- which(
    results$is_dense_support & results$status == "ok"
  )

  dense_row <- if (length(dense_index) > 0L) {
    results[dense_index[1L], , drop = FALSE]
  } else {
    NULL
  }

  selected_fit <- fit_cache[[best_index]]

  if (is.null(selected_fit)) {
    stop("Internal error: selected sparse fit was not retained.", call. = FALSE)
  }

  selected_support <- as.integer(best_row$support_size)

  support_min <- min(results$support_size[results$status == "ok"])
  support_max <- max(results$support_size[results$status == "ok"])

  path_diagnostics <- list(
    selector = selector,
    selected_support = selected_support,
    selected_is_dense = isTRUE(best_row$is_dense_support),
    selected_at_lower_boundary = selected_support == support_min,
    selected_at_upper_boundary = selected_support == support_max,
    support_min = support_min,
    support_max = support_max,
    n_candidates = nrow(results),
    n_successful = sum(results$status == "ok"),
    dense_support_included = any(results$is_dense_support),
    dense_selector_value = if (!is.null(dense_row)) dense_row[[selector]] else NA_real_,
    best_sparse_selector_value = if (!is.null(best_sparse_row)) {
      best_sparse_row[[selector]]
    } else {
      NA_real_
    },
    delta_best_sparse_minus_dense = if (
      !is.null(best_sparse_row) && !is.null(dense_row)
    ) {
      best_sparse_row[[selector]] - dense_row[[selector]]
    } else {
      NA_real_
    }
  )

  out <- list(
    path = results,
    selector = selector,
    selected_row = best_row,
    selected_fit = selected_fit,
    best_sparse_row = best_sparse_row,
    dense_row = dense_row,
    diagnostics = path_diagnostics,
    support_grid = support_grid,
    d = d,
    p = p,
    n = n,
    covariance_method = basis_fit$covariance_method,
    backend = basis_fit$backend,
    dense_basis_fit = basis_fit,
    dense_reference = dense_reference,
    fits = if (isTRUE(store_fits)) fit_cache else NULL,
    call = match.call()
  )

  class(out) <- "ssurvdr_sparse_path"

  out
}
