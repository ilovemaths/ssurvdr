#' Fit reduced Cox proportional hazards model
#'
#' Fits a Cox proportional hazards model using sufficient dimension reduction
#' scores obtained from SDR basis vectors.
#'
#' @param X Numeric predictor matrix with observations in rows.
#' @param time Numeric vector of observed survival times.
#' @param status Numeric or integer censoring indicator vector.
#' @param basis_fit Object returned by [extract_sdr_basis()].
#' @param ties Character string specifying ties handling method for coxph.
#' @param x Logical. Passed to survival::coxph().
#' @param y Logical. Passed to survival::coxph().
#' @param model Logical. Passed to survival::coxph().
#'
#' @return A list containing the reduced predictors, fitted Cox model,
#'   coefficients, linear predictors, likelihood statistics, and diagnostics.
#'
#' @export
fit_reduced_cox <- function(
    X,
    time,
    status,
    basis_fit,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
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

  if (missing(basis_fit)) {
    stop("'basis_fit' must be supplied.", call. = FALSE)
  }

  if (!inherits(basis_fit, "ssurvdr_basis")) {
    stop(
      "'basis_fit' must be an object returned by extract_sdr_basis().",
      call. = FALSE
    )
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

  n <- nrow(X)

  if (length(time) != n) {
    stop("'time' length must equal nrow(X).", call. = FALSE)
  }

  if (length(status) != n) {
    stop("'status' length must equal nrow(X).", call. = FALSE)
  }

  if (anyNA(X)) {
    stop("'X' must not contain missing values.", call. = FALSE)
  }

  if (anyNA(time)) {
    stop("'time' must not contain missing values.", call. = FALSE)
  }

  if (anyNA(status)) {
    stop("'status' must not contain missing values.", call. = FALSE)
  }

  if (any(time <= 0)) {
    stop("'time' must contain strictly positive values.", call. = FALSE)
  }

  status <- as.integer(status)

  if (!all(status %in% c(0L, 1L))) {
    stop("'status' must contain only 0/1 values.", call. = FALSE)
  }

  B <- basis_fit$basis

  if (ncol(X) != nrow(B)) {
    stop(
      "Number of predictors in X must equal number of rows in basis.",
      call. = FALSE
    )
  }

  X_work <- .apply_ssurvdr_preprocessing(
    X = X,
    preprocessing = basis_fit$preprocessing
  )

  Z <- X_work %*% B

  Z <- as.matrix(Z)

  colnames(Z) <- paste0(
    "SDR",
    seq_len(ncol(Z))
  )

  analysis_data <- data.frame(
    time = time,
    status = status,
    Z,
    stringsAsFactors = FALSE
  )

  formula_terms <- paste(
    colnames(Z),
    collapse = " + "
  )

  formula_text <- paste0(
    "survival::Surv(time, status) ~ ",
    formula_terms
  )

  cox_formula <- stats::as.formula(formula_text)

  cox_fit <- survival::coxph(
    formula = cox_formula,
    data = analysis_data,
    ties = ties,
    x = x,
    y = y,
    model = model
  )

  fit_summary <- summary(cox_fit)

  linear_predictor <- stats::predict(
    cox_fit,
    type = "lp"
  )

  risk_score <- stats::predict(
    cox_fit,
    type = "risk"
  )

  concordance <- fit_summary$concordance

  out <- list(
    reduced_predictors = Z,
    model_frame = analysis_data,
    cox_model = cox_fit,
    coefficients = stats::coef(cox_fit),
    hazard_ratios = exp(stats::coef(cox_fit)),
    linear_predictor = linear_predictor,
    risk_score = risk_score,
    loglik = cox_fit$loglik,
    loglik_final = cox_fit$loglik[2],
    aic = stats::AIC(cox_fit),
    n = n,
    d = ncol(Z),
    concordance = concordance,
    basis_fit = basis_fit,
    preprocessing = basis_fit$preprocessing,
    call = match.call()
  )

  class(out) <- "ssurvdr_cox"

  out
}
