#' Fit sparse survival dimension reduction model
#'
#' Fits the current end-to-end ssurvdr workflow: survival-informed slicing,
#' covariance-stabilised survival SAVE kernel construction, SDR basis extraction,
#' reduced Cox proportional hazards modelling, and information-criterion scoring.
#'
#' @param X Numeric predictor matrix with observations in rows.
#' @param time Numeric vector of observed survival times.
#' @param status Numeric or integer censoring indicator vector.
#' @param d Structural dimension.
#' @param n_slices Number of survival-informed slices within each censoring group.
#' @param min_slice_size Minimum number of observations allowed in a slice.
#' @param covariance Character string specifying covariance estimator.
#'   Currently supports "sample" and "mec".
#' @param rho MEC convex mixing weight.
#' @param alpha MEC diagonal stabilisation weight.
#' @param eps Small positive stabilisation constant.
#' @param centre Logical. If TRUE, centre predictors in SAVE construction.
#' @param scale Logical. If TRUE, scale predictors in SAVE construction.
#' @param ties Character string specifying ties handling method for Cox model.
#' @param backend Computational backend passed to [save_kernel()].
#'
#' @return A list containing slicing output, SAVE kernel, SDR basis, reduced Cox
#'   model, information criteria, and diagnostics.
#'
#' @examples
#' dat <- simulate_ssurvdr(n = 100, p = 12, active = 1:3, seed = 2026)
#' fit <- fit_ssurvdr(dat$X, dat$time, dat$status, d = 2,
#'                   n_slices = 3, covariance = "mec", backend = "primal")
#' fit$criteria
#'
#' @export
fit_ssurvdr <- function(
    X,
    time,
    status,
    d = 1,
    n_slices = 5,
    min_slice_size = 2,
    covariance = c("mec", "sample"),
    rho = 0.5,
    alpha = 0.05,
    eps = 1e-8,
    centre = TRUE,
    scale = FALSE,
    ties = "efron",
    backend = c("auto", "primal", "dual")
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

  covariance <- match.arg(covariance)
  backend <- match.arg(backend)

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

  if (length(time) != nrow(X)) {
    stop("'time' length must equal nrow(X).", call. = FALSE)
  }

  if (length(status) != nrow(X)) {
    stop("'status' length must equal nrow(X).", call. = FALSE)
  }

  if (!is.numeric(d) || length(d) != 1L || d < 1) {
    stop("'d' must be a positive integer.", call. = FALSE)
  }

  d <- as.integer(d)

  if (d > ncol(X)) {
    stop("'d' cannot exceed the number of predictors.", call. = FALSE)
  }

  start_time <- Sys.time()

  slices <- survival_slicing(
    time = time,
    status = status,
    n_slices = n_slices,
    min_slice_size = min_slice_size
  )

  save_fit <- save_kernel(
    X = X,
    slices = slices,
    covariance = covariance,
    rho = rho,
    alpha = alpha,
    eps = eps,
    centre = centre,
    scale = scale,
    backend = backend
  )

  basis_fit <- extract_sdr_basis(
    save_fit = save_fit,
    d = d
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

  end_time <- Sys.time()

  runtime_seconds <- as.numeric(
    difftime(
      end_time,
      start_time,
      units = "secs"
    )
  )

  out <- list(
    slices = slices,
    save_fit = save_fit,
    basis_fit = basis_fit,
    cox_fit = cox_fit,
    icomp_fit = icomp_fit,
    criteria = icomp_fit$criteria,
    selected = list(
      covariance = covariance,
      d = d,
      n_slices = n_slices,
      rho = rho,
      alpha = alpha,
      backend = save_fit$backend
    ),
    diagnostics = list(
      n = nrow(X),
      p = ncol(X),
      runtime_seconds = runtime_seconds,
      covariance_method = covariance,
      backend = save_fit$backend,
      row_space_rank = save_fit$row_space_rank,
      data_span_restricted = isTRUE(save_fit$data_span_restricted),
      condition_number = if (!is.null(save_fit$mec_fit)) {
        save_fit$mec_fit$condition_number
      } else {
        NA_real_
      },
      c_index = cox_fit$concordance
    ),
    call = match.call()
  )

  class(out) <- "ssurvdr_fit"

  out
}
