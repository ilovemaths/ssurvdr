#' Summarise ssurvdr model output
#'
#' Produces clean summary tables from objects returned by fit_ssurvdr(),
#' fit_ssurvdr_grid(), and select_variables().
#'
#' @param object An object returned by fit_ssurvdr(), fit_ssurvdr_grid(), or
#'   select_variables().
#' @param top Number of top variables to return when summarising a variable
#'   selection object.
#'
#' @return A list of summary tables.
#'
#' @export
summarise_ssurvdr <- function(
    object,
    top = 10
) {
  if (missing(object)) {
    stop("'object' must be supplied.", call. = FALSE)
  }

  if (!is.numeric(top) || length(top) != 1L || top < 1) {
    stop("'top' must be a positive integer.", call. = FALSE)
  }

  top <- as.integer(top)

  if (inherits(object, "ssurvdr_fit")) {
    model_summary <- data.frame(
      n = object$diagnostics$n,
      p = object$diagnostics$p,
      covariance = object$selected$covariance,
      d = object$selected$d,
      n_slices = object$selected$n_slices,
      rho = object$selected$rho,
      alpha = object$selected$alpha,
      condition_number = object$diagnostics$condition_number,
      runtime_seconds = object$diagnostics$runtime_seconds,
      c_index = as.numeric(object$diagnostics$c_index[1]),
      c_index_se = as.numeric(object$diagnostics$c_index[2]),
      stringsAsFactors = FALSE
    )

    out <- list(
      model_summary = model_summary,
      criteria = object$criteria,
      slice_table = object$slices$slice_table,
      save_diagnostics = object$save_fit$diagnostics
    )

    class(out) <- "ssurvdr_summary"

    return(out)
  }

  if (inherits(object, "ssurvdr_grid")) {
    best <- object$best_model

    best_summary <- data.frame(
      selector = object$selector,
      best_model_id = object$best_model_id,
      covariance = best$selected$covariance,
      d = best$selected$d,
      n_slices = best$selected$n_slices,
      rho = best$selected$rho,
      alpha = best$selected$alpha,
      c_index = as.numeric(best$diagnostics$c_index[1]),
      c_index_se = as.numeric(best$diagnostics$c_index[2]),
      condition_number = best$diagnostics$condition_number,
      runtime_seconds = best$diagnostics$runtime_seconds,
      stringsAsFactors = FALSE
    )

    out <- list(
      best_summary = best_summary,
      ranked_results = object$ranked_results,
      all_results = object$results,
      diagnostics = as.data.frame(
        object$diagnostics,
        stringsAsFactors = FALSE
      )
    )

    class(out) <- "ssurvdr_summary"

    return(out)
  }

  if (inherits(object, "ssurvdr_variable_selection")) {
    n_top <- min(top, nrow(object$scores))

    top_variables <- object$scores[
      seq_len(n_top),
      ,
      drop = FALSE
    ]

    selection_summary <- data.frame(
      p = object$p,
      n_selected = object$n_selected,
      threshold = object$threshold,
      top = ifelse(is.null(object$top), NA_integer_, object$top),
      stringsAsFactors = FALSE
    )

    out <- list(
      selection_summary = selection_summary,
      top_variables = top_variables,
      selected = object$selected
    )

    class(out) <- "ssurvdr_summary"

    return(out)
  }

  stop(
    "'object' must be an ssurvdr_fit, ssurvdr_grid, or ssurvdr_variable_selection object.",
    call. = FALSE
  )
}
