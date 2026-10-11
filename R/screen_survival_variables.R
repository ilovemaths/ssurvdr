#' Screen variables using univariate Cox models
#'
#' Performs preliminary survival variable screening using univariate Cox
#' proportional hazards models. This function is intended as a high-dimensional
#' screening layer before survival SAVE, MEC covariance estimation, and sparse
#' SDR modelling.
#'
#' @param X Numeric predictor matrix with observations in rows.
#' @param time Numeric vector of observed survival times.
#' @param status Numeric or integer censoring indicator vector.
#' @param top Optional number of top variables to retain.
#' @param p_value Optional p-value threshold for screening.
#' @param min_events Minimum number of observed events required.
#' @param ties Ties handling method passed to survival::coxph().
#'
#' @return A list containing screening scores, selected variables, selected
#'   indices, and diagnostics.
#'
#' @export
screen_survival_variables <- function(
    X,
    time,
    status,
    top = NULL,
    p_value = NULL,
    min_events = 5,
    ties = "efron"
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

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    stop("'X' must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (!is.numeric(X)) {
    stop("'X' must contain numeric values only.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' must contain only finite non-missing values.", call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)

  if (length(time) != n) {
    stop("'time' length must equal nrow(X).", call. = FALSE)
  }

  if (length(status) != n) {
    stop("'status' length must equal nrow(X).", call. = FALSE)
  }

  if (anyNA(time) || any(!is.finite(time))) {
    stop("'time' must contain only finite non-missing values.", call. = FALSE)
  }

  if (any(time <= 0)) {
    stop("'time' must contain strictly positive values.", call. = FALSE)
  }

  status <- as.integer(status)

  if (anyNA(status) || !all(status %in% c(0L, 1L))) {
    stop("'status' must contain only 0/1 values.", call. = FALSE)
  }

  if (sum(status == 1L) < min_events) {
    stop("The number of observed events is smaller than 'min_events'.", call. = FALSE)
  }

  if (!is.null(top)) {
    if (!is.numeric(top) || length(top) != 1L || top < 1) {
      stop("'top' must be NULL or a positive integer.", call. = FALSE)
    }

    top <- as.integer(top)

    if (top > p) {
      stop("'top' cannot exceed the number of predictors.", call. = FALSE)
    }
  }

  if (!is.null(p_value)) {
    if (!is.numeric(p_value) || length(p_value) != 1L || p_value <= 0 || p_value > 1) {
      stop("'p_value' must be NULL or a numeric value in (0, 1].", call. = FALSE)
    }
  }

  variable_names <- colnames(X)

  if (is.null(variable_names)) {
    variable_names <- paste0("X", seq_len(p))
  }

  results <- data.frame(
    variable = variable_names,
    index = seq_len(p),
    coefficient = NA_real_,
    hazard_ratio = NA_real_,
    se = NA_real_,
    z = NA_real_,
    p_value = NA_real_,
    score = NA_real_,
    converged = FALSE,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(p)) {
    xj <- X[, j]

    if (stats::sd(xj) <= .Machine$double.eps) {
      results$error_message[j] <- "near-zero variance"
      next
    }

    dat <- data.frame(
      time = time,
      status = status,
      x = xj
    )

    fit_j <- tryCatch(
      {
        survival::coxph(
          survival::Surv(time, status) ~ x,
          data = dat,
          ties = ties
        )
      },
      error = function(e) {
        e
      },
      warning = function(w) {
        suppressWarnings(
          survival::coxph(
            survival::Surv(time, status) ~ x,
            data = dat,
            ties = ties
          )
        )
      }
    )

    if (inherits(fit_j, "error")) {
      results$error_message[j] <- conditionMessage(fit_j)
      next
    }

    summary_j <- summary(fit_j)

    coef_table <- summary_j$coefficients

    results$coefficient[j] <- coef_table[1, "coef"]
    results$hazard_ratio[j] <- coef_table[1, "exp(coef)"]
    results$se[j] <- coef_table[1, "se(coef)"]
    results$z[j] <- coef_table[1, "z"]
    results$p_value[j] <- coef_table[1, "Pr(>|z|)"]
    results$score[j] <- abs(coef_table[1, "z"])
    results$converged[j] <- TRUE
  }

  results$rank <- rank(
    -results$score,
    ties.method = "first",
    na.last = "keep"
  )

  results <- results[
    order(results$rank),
    ,
    drop = FALSE
  ]

  rownames(results) <- NULL

  valid_results <- results[is.finite(results$score), , drop = FALSE]

  if (nrow(valid_results) < 1L) {
    stop("No variable could be screened successfully.", call. = FALSE)
  }

  if (!is.null(top)) {
    selected <- valid_results[
      seq_len(min(top, nrow(valid_results))),
      ,
      drop = FALSE
    ]
  } else if (!is.null(p_value)) {
    selected <- valid_results[
      valid_results$p_value <= p_value,
      ,
      drop = FALSE
    ]
  } else {
    selected <- valid_results
  }

  out <- list(
    scores = results,
    selected = selected,
    selected_variables = selected$variable,
    selected_indices = selected$index,
    top = top,
    p_value = p_value,
    n_selected = nrow(selected),
    n_screened = p,
    n_successful = sum(results$converged),
    n_failed = sum(!results$converged),
    min_events = min_events
  )

  class(out) <- "ssurvdr_screening"

  out
}
