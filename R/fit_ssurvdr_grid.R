#' Fit ssurvdr model over a tuning grid
#'
#' Fits the ssurvdr workflow over candidate covariance estimators and structural
#' dimensions, then returns model-comparison results using AIC, BIC, CAIC,
#' ICOMP(IFIM), ICOMP(misspec), C-index, runtime, and condition number.
#'
#' The survival slices and SAVE operator are computed once per covariance
#' estimator and reused across structural dimensions, avoiding redundant kernel
#' construction and eigendecomposition.
#'
#' @param X Numeric predictor matrix with observations in rows.
#' @param time Numeric vector of observed survival times.
#' @param status Numeric or integer censoring indicator vector.
#' @param d_grid Integer vector of candidate structural dimensions.
#' @param covariance_grid Character vector of covariance estimators.
#'   Currently supports "mec" and "sample".
#' @param n_slices Number of survival-informed slices within each censoring group.
#' @param min_slice_size Minimum number of observations allowed in a slice.
#' @param rho MEC convex mixing weight.
#' @param alpha MEC diagonal stabilisation weight.
#' @param eps Small positive stabilisation constant.
#' @param centre Logical. If TRUE, centre predictors in SAVE construction.
#' @param scale Logical. If TRUE, scale predictors in SAVE construction.
#' @param ties Character string specifying ties handling method for Cox model.
#' @param selector Criterion used to rank models. One of "ICOMP_IFIM",
#'   "ICOMP_MISSPEC", "AIC", "BIC", or "CAIC".
#' @param continue_on_error Logical. If TRUE, failed candidate models are
#'   recorded and the grid search continues.
#' @param backend Computational backend passed to [save_kernel()].
#' @param reuse_kernel Logical. If TRUE, compute the SAVE operator once per
#'   covariance estimator and reuse it across d. FALSE retains the legacy
#'   fit-by-fit behaviour and is mainly intended for regression testing.
#'
#' @return A list containing all fitted models, a model-comparison table, the
#'   best model, cached SAVE operators, and diagnostics.
#'
#' @export
fit_ssurvdr_grid <- function(
    X,
    time,
    status,
    d_grid = 1:3,
    covariance_grid = c("mec", "sample"),
    n_slices = 5,
    min_slice_size = 2,
    rho = 0.5,
    alpha = 0.05,
    eps = 1e-8,
    centre = TRUE,
    scale = FALSE,
    ties = "efron",
    selector = c(
      "ICOMP_IFIM",
      "ICOMP_MISSPEC",
      "AIC",
      "BIC",
      "CAIC"
    ),
    continue_on_error = TRUE,
    backend = c("auto", "primal", "dual"),
    reuse_kernel = TRUE
) {
  if (missing(X)) {
    stop("'X' must be supplied.", call. = FALSE)
  }

  if (missing(time)) {
    stop("'time' must be supplied.", call. = FALSE)
  }

  if (missing(status)) {
    stop("'status' must be supplied.", call. = FALSE)
  }

  selector <- match.arg(selector)
  backend <- match.arg(backend)

  if (!is.logical(reuse_kernel) || length(reuse_kernel) != 1L) {
    stop("'reuse_kernel' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("'X' must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' must contain only finite non-missing values.", call. = FALSE)
  }

  if (length(time) != nrow(X)) {
    stop("'time' length must equal nrow(X).", call. = FALSE)
  }

  if (length(status) != nrow(X)) {
    stop("'status' length must equal nrow(X).", call. = FALSE)
  }

  allowed_covariance <- c("mec", "sample")

  if (!all(covariance_grid %in% allowed_covariance)) {
    stop(
      "'covariance_grid' currently supports only 'mec' and 'sample'.",
      call. = FALSE
    )
  }

  if (!is.numeric(d_grid) || length(d_grid) < 1L) {
    stop("'d_grid' must be a non-empty numeric or integer vector.", call. = FALSE)
  }

  d_grid <- sort(unique(as.integer(d_grid)))

  if (any(d_grid < 1L)) {
    stop("All values in 'd_grid' must be positive integers.", call. = FALSE)
  }

  if (max(d_grid) > ncol(X)) {
    stop("No value in 'd_grid' may exceed the number of predictors.", call. = FALSE)
  }

  grid <- expand.grid(
    covariance = covariance_grid,
    d = d_grid,
    stringsAsFactors = FALSE
  )

  fits <- vector("list", nrow(grid))
  kernel_cache <- vector("list", length(covariance_grid))
  names(kernel_cache) <- covariance_grid
  kernel_runtime <- stats::setNames(rep(NA_real_, length(covariance_grid)), covariance_grid)
  kernel_error <- stats::setNames(rep(NA_character_, length(covariance_grid)), covariance_grid)

  results <- data.frame(
    model_id = seq_len(nrow(grid)),
    covariance = grid$covariance,
    d = grid$d,
    AIC = NA_real_,
    BIC = NA_real_,
    CAIC = NA_real_,
    ICOMP_IFIM = NA_real_,
    ICOMP_MISSPEC = NA_real_,
    c_index = NA_real_,
    c_index_se = NA_real_,
    condition_number = NA_real_,
    runtime_seconds = NA_real_,
    kernel_runtime_seconds = NA_real_,
    downstream_runtime_seconds = NA_real_,
    backend = NA_character_,
    row_space_rank = NA_integer_,
    status = "not_run",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  overall_start <- Sys.time()

  # Legacy path retained for regression checks.
  if (!isTRUE(reuse_kernel)) {
    for (i in seq_len(nrow(grid))) {
      start_i <- Sys.time()

      fit_i <- tryCatch(
        fit_ssurvdr(
          X = X,
          time = time,
          status = status,
          d = grid$d[i],
          n_slices = n_slices,
          min_slice_size = min_slice_size,
          covariance = grid$covariance[i],
          rho = rho,
          alpha = alpha,
          eps = eps,
          centre = centre,
          scale = scale,
          ties = ties,
          backend = backend
        ),
        error = function(e) e
      )

      if (inherits(fit_i, "error")) {
        results$status[i] <- "failed"
        results$error_message[i] <- conditionMessage(fit_i)

        if (!isTRUE(continue_on_error)) {
          stop(results$error_message[i], call. = FALSE)
        }
        next
      }

      fits[[i]] <- fit_i
      crit <- fit_i$criteria

      results$AIC[i] <- crit$value[crit$criterion == "AIC"]
      results$BIC[i] <- crit$value[crit$criterion == "BIC"]
      results$CAIC[i] <- crit$value[crit$criterion == "CAIC"]
      results$ICOMP_IFIM[i] <- crit$value[crit$criterion == "ICOMP_IFIM"]
      results$ICOMP_MISSPEC[i] <- crit$value[crit$criterion == "ICOMP_MISSPEC"]

      concordance <- fit_i$cox_fit$concordance
      if (length(concordance) >= 1L) results$c_index[i] <- as.numeric(concordance[1])
      if (length(concordance) >= 2L) results$c_index_se[i] <- as.numeric(concordance[2])

      results$condition_number[i] <- fit_i$diagnostics$condition_number
      results$runtime_seconds[i] <- as.numeric(difftime(Sys.time(), start_i, units = "secs"))
      results$backend[i] <- fit_i$diagnostics$backend
      results$row_space_rank[i] <- fit_i$diagnostics$row_space_rank
      results$status[i] <- "ok"
    }
  } else {
    slices <- survival_slicing(
      time = time,
      status = status,
      n_slices = n_slices,
      min_slice_size = min_slice_size
    )

    # Compute each SAVE operator exactly once.
    for (covariance_method in covariance_grid) {
      kernel_start <- Sys.time()

      kernel_i <- tryCatch(
        save_kernel(
          X = X,
          slices = slices,
          covariance = covariance_method,
          rho = rho,
          alpha = alpha,
          eps = eps,
          centre = centre,
          scale = scale,
          backend = backend
        ),
        error = function(e) e
      )

      kernel_runtime[covariance_method] <- as.numeric(
        difftime(Sys.time(), kernel_start, units = "secs")
      )

      if (inherits(kernel_i, "error")) {
        kernel_error[covariance_method] <- conditionMessage(kernel_i)
      } else {
        kernel_cache[[covariance_method]] <- kernel_i
      }
    }

    for (i in seq_len(nrow(grid))) {
      current_covariance <- grid$covariance[i]
      current_d <- grid$d[i]
      save_fit <- kernel_cache[[current_covariance]]

      results$kernel_runtime_seconds[i] <- kernel_runtime[current_covariance]

      if (is.null(save_fit)) {
        results$status[i] <- "failed"
        results$error_message[i] <- kernel_error[current_covariance]

        if (!isTRUE(continue_on_error)) {
          stop(results$error_message[i], call. = FALSE)
        }
        next
      }

      downstream_start <- Sys.time()

      fit_i <- tryCatch(
        {
          basis_fit <- extract_sdr_basis(
            save_fit = save_fit,
            d = current_d
          )

          cox_fit <- fit_reduced_cox(
            X = X,
            time = time,
            status = status,
            basis_fit = basis_fit,
            ties = ties
          )

          icomp_fit <- compute_icomp(
            cox_fit = cox_fit,
            eps = eps
          )

          downstream_seconds <- as.numeric(
            difftime(Sys.time(), downstream_start, units = "secs")
          )

          out <- list(
            slices = slices,
            save_fit = save_fit,
            basis_fit = basis_fit,
            cox_fit = cox_fit,
            icomp_fit = icomp_fit,
            criteria = icomp_fit$criteria,
            selected = list(
              covariance = current_covariance,
              d = current_d,
              n_slices = n_slices,
              rho = rho,
              alpha = alpha,
              backend = save_fit$backend
            ),
            diagnostics = list(
              n = nrow(X),
              p = ncol(X),
              runtime_seconds = kernel_runtime[current_covariance] + downstream_seconds,
              kernel_runtime_seconds = kernel_runtime[current_covariance],
              downstream_runtime_seconds = downstream_seconds,
              covariance_method = current_covariance,
              backend = save_fit$backend,
              row_space_rank = save_fit$row_space_rank,
              data_span_restricted = isTRUE(save_fit$data_span_restricted),
              condition_number = if (!is.null(save_fit$mec_fit)) {
                save_fit$mec_fit$condition_number
              } else {
                NA_real_
              },
              c_index = cox_fit$concordance,
              kernel_reused = TRUE
            ),
            call = match.call()
          )

          class(out) <- "ssurvdr_fit"
          out
        },
        error = function(e) e
      )

      if (inherits(fit_i, "error")) {
        results$status[i] <- "failed"
        results$error_message[i] <- conditionMessage(fit_i)

        if (!isTRUE(continue_on_error)) {
          stop(results$error_message[i], call. = FALSE)
        }
        next
      }

      fits[[i]] <- fit_i
      crit <- fit_i$criteria

      results$AIC[i] <- crit$value[crit$criterion == "AIC"]
      results$BIC[i] <- crit$value[crit$criterion == "BIC"]
      results$CAIC[i] <- crit$value[crit$criterion == "CAIC"]
      results$ICOMP_IFIM[i] <- crit$value[crit$criterion == "ICOMP_IFIM"]
      results$ICOMP_MISSPEC[i] <- crit$value[crit$criterion == "ICOMP_MISSPEC"]

      concordance <- fit_i$cox_fit$concordance
      if (length(concordance) >= 1L) results$c_index[i] <- as.numeric(concordance[1])
      if (length(concordance) >= 2L) results$c_index_se[i] <- as.numeric(concordance[2])

      results$condition_number[i] <- fit_i$diagnostics$condition_number
      results$runtime_seconds[i] <- fit_i$diagnostics$runtime_seconds
      results$downstream_runtime_seconds[i] <- fit_i$diagnostics$downstream_runtime_seconds
      results$backend[i] <- fit_i$diagnostics$backend
      results$row_space_rank[i] <- fit_i$diagnostics$row_space_rank
      results$status[i] <- "ok"
    }
  }

  valid_results <- results[results$status == "ok", , drop = FALSE]

  if (nrow(valid_results) < 1L) {
    stop("No model in the grid search completed successfully.", call. = FALSE)
  }

  ranked_results <- valid_results[
    order(valid_results[[selector]], decreasing = FALSE),
    ,
    drop = FALSE
  ]
  rownames(ranked_results) <- NULL

  best_model_id <- ranked_results$model_id[1]
  best_fit <- fits[[best_model_id]]

  total_runtime_seconds <- as.numeric(
    difftime(Sys.time(), overall_start, units = "secs")
  )

  out <- list(
    results = results,
    ranked_results = ranked_results,
    best_model = best_fit,
    best_model_id = best_model_id,
    selector = selector,
    grid = grid,
    fits = fits,
    kernel_cache = kernel_cache,
    diagnostics = list(
      n_models = nrow(grid),
      n_successful = sum(results$status == "ok"),
      n_failed = sum(results$status == "failed"),
      total_runtime_seconds = total_runtime_seconds,
      backend_requested = backend,
      reuse_kernel = isTRUE(reuse_kernel),
      kernel_runtime_seconds = kernel_runtime,
      kernel_error = kernel_error
    ),
    call = match.call()
  )

  class(out) <- "ssurvdr_grid"
  out
}
