#' Select variables from SDR basis
#'
#' Selects variables using row norms of the estimated sufficient dimension
#' reduction basis matrix.
#'
#' @param basis_fit Object returned by [extract_sdr_basis()] or a numeric basis matrix.
#' @param threshold Numeric threshold for variable selection.
#' @param top Optional integer. If supplied, returns the top variables by row norm.
#' @param variable_names Optional character vector of predictor names.
#'
#' @return A list containing variable scores, selected variables, selected indices,
#'   and selection diagnostics.
#'
#' @export
select_variables <- function(
    basis_fit,
    threshold = 1e-6,
    top = NULL,
    variable_names = NULL
) {
  if (missing(basis_fit)) {
    stop("'basis_fit' must be supplied.", call. = FALSE)
  }

  if (inherits(basis_fit, "ssurvdr_basis")) {
    B <- basis_fit$basis
  } else {
    B <- basis_fit
  }

  if (!is.matrix(B)) {
    stop("'basis_fit' must be an ssurvdr_basis object or a numeric matrix.", call. = FALSE)
  }

  if (!is.numeric(B)) {
    stop("The SDR basis matrix must be numeric.", call. = FALSE)
  }

  if (anyNA(B)) {
    stop("The SDR basis matrix must not contain missing values.", call. = FALSE)
  }

  p <- nrow(B)

  if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0) {
    stop("'threshold' must be a single non-negative numeric value.", call. = FALSE)
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

  if (is.null(variable_names)) {
    variable_names <- rownames(B)

    if (is.null(variable_names)) {
      variable_names <- paste0("V", seq_len(p))
    }
  }

  if (length(variable_names) != p) {
    stop("'variable_names' must have length equal to the number of predictors.", call. = FALSE)
  }

  row_norms <- sqrt(rowSums(B^2))

  ranking <- order(row_norms, decreasing = TRUE)

  scores <- data.frame(
    variable = variable_names,
    index = seq_len(p),
    score = row_norms,
    rank = rank(
      -row_norms,
      ties.method = "first"
    ),
    selected_by_threshold = row_norms > threshold,
    stringsAsFactors = FALSE
  )

  scores <- scores[order(scores$rank), , drop = FALSE]
  rownames(scores) <- NULL

  if (is.null(top)) {
    selected_scores <- scores[scores$selected_by_threshold, , drop = FALSE]
  } else {
    selected_scores <- scores[seq_len(top), , drop = FALSE]
  }

  out <- list(
    scores = scores,
    selected = selected_scores,
    selected_variables = selected_scores$variable,
    selected_indices = selected_scores$index,
    threshold = threshold,
    top = top,
    n_selected = nrow(selected_scores),
    p = p
  )

  class(out) <- "ssurvdr_variable_selection"

  out
}
