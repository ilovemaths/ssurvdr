#' Evaluate ssurvdr fitted model against known simulation truth
#'
#' Evaluates selected variables, subspace recovery, and survival prediction
#' performance when the true active set and true SDR basis are known.
#'
#' @param fit Object returned by [fit_ssurvdr()] or [fit_ssurvdr_grid()].
#' @param simulation Object returned by [simulate_ssurvdr()].
#' @param variable_selection Optional object returned by [select_variables()].
#'   If NULL, variable selection is computed from the fitted basis.
#' @param top Optional number of variables to select by row-norm ranking.
#' @param threshold Numeric threshold for row-norm selection when top is NULL.
#' @param screen_indices Optional integer vector giving the original predictor
#'   indices retained after screening. If supplied, the estimated basis is
#'   expanded back to the original predictor dimension before projection
#'   distance is computed. If `variable_selection` is supplied, its selected
#'   indices must already refer to the original predictor matrix.
#'
#' @return A list containing variable-selection metrics, projection distance,
#'   C-index, selected variables, and diagnostics.
#'
#' @export
evaluate_ssurvdr <- function(
    fit,
    simulation,
    variable_selection = NULL,
    top = NULL,
    threshold = 1e-6,
    screen_indices = NULL
) {
  if (missing(fit)) {
    stop("'fit' must be supplied.", call. = FALSE)
  }

  if (missing(simulation)) {
    stop("'simulation' must be supplied.", call. = FALSE)
  }

  if (!inherits(simulation, "ssurvdr_simulation")) {
    stop(
      "'simulation' must be an object returned by simulate_ssurvdr().",
      call. = FALSE
    )
  }

  if (inherits(fit, "ssurvdr_grid")) {
    fit_object <- fit$best_model
  } else if (inherits(fit, "ssurvdr_fit")) {
    fit_object <- fit
  } else {
    stop(
      "'fit' must be an ssurvdr_fit or ssurvdr_grid object.",
      call. = FALSE
    )
  }

  if (!is.null(screen_indices)) {
    if (!is.numeric(screen_indices) || anyNA(screen_indices) ||
        any(!is.finite(screen_indices)) ||
        any(screen_indices != as.integer(screen_indices)) ||
        anyDuplicated(screen_indices) ||
        any(screen_indices < 1L | screen_indices > ncol(simulation$X)) ||
        length(screen_indices) != nrow(fit_object$basis_fit$basis)) {
      stop("'screen_indices' must uniquely map each fitted predictor to an original predictor.",
           call. = FALSE)
    }
    screen_indices <- as.integer(screen_indices)
  }

  if (is.null(variable_selection)) {
    variable_names <- colnames(simulation$X)
    if (!is.null(screen_indices) && !is.null(variable_names)) {
      variable_names <- variable_names[screen_indices]
    }
    variable_selection <- select_variables(
      basis_fit = fit_object$basis_fit,
      threshold = threshold,
      top = top,
      variable_names = variable_names
    )
    if (!is.null(screen_indices)) {
      original_indices <- screen_indices[variable_selection$selected_indices]
      variable_selection$selected_indices <- original_indices
      variable_selection$selected$index <- original_indices
    }
  }

  if (!inherits(variable_selection, "ssurvdr_variable_selection")) {
    stop(
      "'variable_selection' must be NULL or an object returned by select_variables().",
      call. = FALSE
    )
  }

  true_active <- sort(unique(as.integer(simulation$active)))
  selected_active <- sort(unique(as.integer(variable_selection$selected_indices)))

  p <- ncol(simulation$X)

  true_inactive <- setdiff(seq_len(p), true_active)
  selected_inactive <- setdiff(seq_len(p), selected_active)

  true_positive <- length(intersect(selected_active, true_active))
  false_positive <- length(setdiff(selected_active, true_active))
  false_negative <- length(setdiff(true_active, selected_active))
  true_negative <- length(intersect(selected_inactive, true_inactive))

  sensitivity <- if (length(true_active) > 0L) {
    true_positive / length(true_active)
  } else {
    NA_real_
  }

  specificity <- if (length(true_inactive) > 0L) {
    true_negative / length(true_inactive)
  } else {
    NA_real_
  }

  precision <- if (length(selected_active) > 0L) {
    true_positive / length(selected_active)
  } else {
    NA_real_
  }

  false_discovery_rate <- if (length(selected_active) > 0L) {
    false_positive / length(selected_active)
  } else {
    NA_real_
  }

  estimated_basis <- fit_object$basis_fit$basis

  if (!is.null(screen_indices)) {
    screen_indices <- as.integer(screen_indices)

    if (any(screen_indices < 1L) || any(screen_indices > ncol(simulation$X))) {
      stop(
        "'screen_indices' must contain valid original predictor indices.",
        call. = FALSE
      )
    }

    if (length(screen_indices) != nrow(estimated_basis)) {
      stop(
        "Length of 'screen_indices' must equal the number of rows in the estimated basis.",
        call. = FALSE
      )
    }

    expanded_basis <- matrix(
      0,
      nrow = ncol(simulation$X),
      ncol = ncol(estimated_basis)
    )

    expanded_basis[screen_indices, ] <- estimated_basis

    estimated_basis <- expanded_basis
  }
  true_basis <- simulation$true_basis

  if (!is.matrix(true_basis)) {
    true_basis <- as.matrix(true_basis)
  }

  projection_matrix <- function(B) {
    B <- as.matrix(B)

    qr_B <- qr(B)

    Q <- qr.Q(qr_B)

    Q %*% t(Q)
  }

  P_hat <- projection_matrix(estimated_basis)
  P_true <- projection_matrix(true_basis)

  projection_distance <- sqrt(sum((P_hat - P_true)^2))

  c_index <- as.numeric(fit_object$cox_fit$concordance[1])
  c_index_se <- as.numeric(fit_object$cox_fit$concordance[2])

  metrics <- data.frame(
    n = nrow(simulation$X),
    p = p,
    n_true_active = length(true_active),
    n_selected = length(selected_active),
    true_positive = true_positive,
    false_positive = false_positive,
    false_negative = false_negative,
    true_negative = true_negative,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    false_discovery_rate = false_discovery_rate,
    projection_distance = projection_distance,
    c_index = c_index,
    c_index_se = c_index_se,
    stringsAsFactors = FALSE
  )

  out <- list(
    metrics = metrics,
    selected_variables = variable_selection$selected,
    true_active_indices = true_active,
    selected_indices = selected_active,
    true_positive_indices = intersect(selected_active, true_active),
    false_positive_indices = setdiff(selected_active, true_active),
    false_negative_indices = setdiff(true_active, selected_active),
    projection_distance = projection_distance,
    fit = fit_object,
    simulation = simulation,
    variable_selection = variable_selection
  )

  class(out) <- "ssurvdr_evaluation"

  out
}
