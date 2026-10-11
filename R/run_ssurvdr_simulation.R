#' Run Monte Carlo simulation for ssurvdr
#'
#' Runs repeated simulation experiments for sparse survival dimension reduction
#' under low-dimensional, borderline-dimensional, high-dimensional, and
#' ultra-high-dimensional settings.
#'
#' @param n Sample size.
#' @param p Number of predictors.
#' @param active Integer vector of active predictor indices.
#' @param covariance Data-generating covariance structure.
#' @param rho Correlation parameter.
#' @param signal_strength Survival signal strength.
#' @param censoring_rate Exponential censoring rate.
#' @param n_rep Number of Monte Carlo replications.
#' @param screening Logical. If TRUE, univariate Cox screening is performed.
#' @param screen_top Number of variables retained after screening.
#' @param covariance_methods Covariance methods used in fit_ssurvdr_grid().
#' @param d_grid Candidate structural dimensions.
#' @param n_slices Number of survival-informed slices.
#' @param selector Model selection criterion.
#' @param variable_top Number of variables selected for evaluation.
#' @param seed Optional random seed.
#' @param continue_on_error Logical. If TRUE, failed replications are recorded.
#'
#' @return A list containing replication-level results, summary results, and
#'   simulation settings.
#'
#' @export
run_ssurvdr_simulation <- function(
    n = 120,
    p = 30,
    active = 1:5,
    covariance = c("identity", "ar1", "compound", "block"),
    rho = 0.5,
    signal_strength = 1.2,
    censoring_rate = 0.08,
    n_rep = 50,
    screening = TRUE,
    screen_top = NULL,
    covariance_methods = c("mec", "sample"),
    d_grid = 1:3,
    n_slices = 4,
    selector = c(
      "ICOMP_IFIM",
      "ICOMP_MISSPEC",
      "AIC",
      "BIC",
      "CAIC"
    ),
    variable_top = NULL,
    seed = NULL,
    continue_on_error = TRUE
) {
  covariance <- match.arg(covariance)
  selector <- match.arg(selector)

  if (!is.numeric(n) || length(n) != 1L || n < 2) {
    stop("'n' must be a single integer greater than or equal to 2.", call. = FALSE)
  }

  if (!is.numeric(p) || length(p) != 1L || p < 1) {
    stop("'p' must be a single positive integer.", call. = FALSE)
  }

  n <- as.integer(n)
  p <- as.integer(p)

  if (!is.numeric(n_rep) || length(n_rep) != 1L || n_rep < 1) {
    stop("'n_rep' must be a positive integer.", call. = FALSE)
  }

  n_rep <- as.integer(n_rep)

  if (!is.numeric(active) || length(active) < 1L) {
    stop("'active' must be a non-empty integer vector.", call. = FALSE)
  }

  active <- sort(unique(as.integer(active)))

  if (any(active < 1L) || any(active > p)) {
    stop("All entries of 'active' must be valid predictor indices.", call. = FALSE)
  }

  if (is.null(screen_top)) {
    screen_top <- min(p, max(ceiling(n / 3), length(active) * 3L))
  }

  if (!is.numeric(screen_top) || length(screen_top) != 1L || screen_top < 1) {
    stop("'screen_top' must be NULL or a positive integer.", call. = FALSE)
  }

  screen_top <- as.integer(screen_top)

  if (screen_top > p) {
    stop("'screen_top' cannot exceed 'p'.", call. = FALSE)
  }

  if (is.null(variable_top)) {
    variable_top <- min(10L, screen_top, p)
  }

  if (!is.numeric(variable_top) || length(variable_top) != 1L || variable_top < 1) {
    stop("'variable_top' must be NULL or a positive integer.", call. = FALSE)
  }

  variable_top <- as.integer(variable_top)

  if (variable_top > p) {
    stop("'variable_top' cannot exceed 'p'.", call. = FALSE)
  }

  dimension_regime <- if (p < n) {
    "low-dimensional"
  } else if (p == n) {
    "borderline-dimensional"
  } else if (p <= 2 * n) {
    "high-dimensional"
  } else {
    "ultra-high-dimensional"
  }

  if (!is.null(seed)) {
    set.seed(seed)
    replication_seeds <- sample.int(.Machine$integer.max, n_rep)
  } else {
    replication_seeds <- sample.int(.Machine$integer.max, n_rep)
  }

  replication_results <- vector("list", n_rep)

  result_table <- data.frame(
    replication = seq_len(n_rep),
    status = "not_run",
    error_message = NA_character_,
    n = n,
    p = p,
    dimension_regime = dimension_regime,
    covariance_data = covariance,
    rho = rho,
    signal_strength = signal_strength,
    censoring_rate = censoring_rate,
    censoring_proportion = NA_real_,
    screening = screening,
    screen_top = screen_top,
    n_screened = NA_integer_,
    best_covariance = NA_character_,
    best_d = NA_integer_,
    selector = selector,
    best_selector_value = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    CAIC = NA_real_,
    ICOMP_IFIM = NA_real_,
    ICOMP_MISSPEC = NA_real_,
    c_index = NA_real_,
    c_index_se = NA_real_,
    condition_number = NA_real_,
    runtime_seconds = NA_real_,
    n_selected = NA_integer_,
    true_positive = NA_integer_,
    false_positive = NA_integer_,
    false_negative = NA_integer_,
    true_negative = NA_integer_,
    sensitivity = NA_real_,
    specificity = NA_real_,
    precision = NA_real_,
    false_discovery_rate = NA_real_,
    projection_distance = NA_real_,
    stringsAsFactors = FALSE
  )

  overall_start <- Sys.time()

  for (b in seq_len(n_rep)) {
    fit_b <- tryCatch(
      {
        sim_b <- simulate_ssurvdr(
          n = n,
          p = p,
          active = active,
          covariance = covariance,
          rho = rho,
          signal_strength = signal_strength,
          censoring_rate = censoring_rate,
          seed = replication_seeds[b]
        )

        if (isTRUE(screening)) {
          screen_b <- screen_survival_variables(
            X = sim_b$X,
            time = sim_b$time,
            status = sim_b$status,
            top = screen_top
          )

          X_model <- sim_b$X[, screen_b$selected_indices, drop = FALSE]
          screen_indices <- screen_b$selected_indices
        } else {
          screen_b <- NULL
          X_model <- sim_b$X
          screen_indices <- NULL
        }

        grid_b <- fit_ssurvdr_grid(
          X = X_model,
          time = sim_b$time,
          status = sim_b$status,
          d_grid = d_grid,
          covariance_grid = covariance_methods,
          n_slices = n_slices,
          selector = selector,
          centre = TRUE,
          scale = TRUE
        )

        var_b <- select_variables(
          basis_fit = grid_b$best_model$basis_fit,
          top = min(variable_top, ncol(X_model)),
          variable_names = colnames(X_model)
        )

        if (isTRUE(screening)) {
          selected_original_indices <- screen_indices[
            match(var_b$selected_variables, colnames(X_model))
          ]

          var_original <- var_b
          var_original$selected$index <- selected_original_indices
          var_original$selected_indices <- selected_original_indices
        } else {
          var_original <- var_b
        }

        eval_b <- evaluate_ssurvdr(
          fit = grid_b,
          simulation = sim_b,
          variable_selection = var_original,
          screen_indices = screen_indices
        )

        list(
          simulation = sim_b,
          screening = screen_b,
          grid = grid_b,
          variable_selection = var_original,
          evaluation = eval_b
        )
      },
      error = function(e) {
        e
      }
    )

    if (inherits(fit_b, "error")) {
      result_table$status[b] <- "failed"
      result_table$error_message[b] <- conditionMessage(fit_b)

      replication_results[[b]] <- NULL

      if (!isTRUE(continue_on_error)) {
        stop(
          paste0(
            "Simulation stopped at replication ",
            b,
            ": ",
            conditionMessage(fit_b)
          ),
          call. = FALSE
        )
      }

      next
    }

    replication_results[[b]] <- fit_b

    best_model <- fit_b$grid$best_model
    criteria <- best_model$criteria
    eval_metrics <- fit_b$evaluation$metrics

    result_table$status[b] <- "ok"
    result_table$censoring_proportion[b] <- fit_b$simulation$censoring_proportion
    result_table$n_screened[b] <- ncol(
      if (isTRUE(screening)) {
        fit_b$simulation$X[, fit_b$screening$selected_indices, drop = FALSE]
      } else {
        fit_b$simulation$X
      }
    )

    result_table$best_covariance[b] <- best_model$selected$covariance
    result_table$best_d[b] <- best_model$selected$d

    result_table$AIC[b] <- criteria$value[criteria$criterion == "AIC"]
    result_table$BIC[b] <- criteria$value[criteria$criterion == "BIC"]
    result_table$CAIC[b] <- criteria$value[criteria$criterion == "CAIC"]
    result_table$ICOMP_IFIM[b] <- criteria$value[criteria$criterion == "ICOMP_IFIM"]
    result_table$ICOMP_MISSPEC[b] <- criteria$value[criteria$criterion == "ICOMP_MISSPEC"]
    result_table$best_selector_value[b] <- result_table[[selector]][b]

    result_table$c_index[b] <- eval_metrics$c_index
    result_table$c_index_se[b] <- eval_metrics$c_index_se
    result_table$condition_number[b] <- best_model$diagnostics$condition_number
    result_table$runtime_seconds[b] <- best_model$diagnostics$runtime_seconds

    result_table$n_selected[b] <- eval_metrics$n_selected
    result_table$true_positive[b] <- eval_metrics$true_positive
    result_table$false_positive[b] <- eval_metrics$false_positive
    result_table$false_negative[b] <- eval_metrics$false_negative
    result_table$true_negative[b] <- eval_metrics$true_negative
    result_table$sensitivity[b] <- eval_metrics$sensitivity
    result_table$specificity[b] <- eval_metrics$specificity
    result_table$precision[b] <- eval_metrics$precision
    result_table$false_discovery_rate[b] <- eval_metrics$false_discovery_rate
    result_table$projection_distance[b] <- eval_metrics$projection_distance
  }

  successful <- result_table[result_table$status == "ok", , drop = FALSE]

  summarise_numeric <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE)
    )
  }

  if (nrow(successful) > 0L) {
    summary_table <- data.frame(
      n = n,
      p = p,
      dimension_regime = dimension_regime,
      covariance_data = covariance,
      rho = rho,
      signal_strength = signal_strength,
      censoring_rate = censoring_rate,
      screening = screening,
      screen_top = screen_top,
      n_rep = n_rep,
      n_successful = nrow(successful),
      success_rate = nrow(successful) / n_rep,
      mean_censoring_proportion = mean(successful$censoring_proportion, na.rm = TRUE),
      mean_sensitivity = mean(successful$sensitivity, na.rm = TRUE),
      sd_sensitivity = stats::sd(successful$sensitivity, na.rm = TRUE),
      mean_specificity = mean(successful$specificity, na.rm = TRUE),
      sd_specificity = stats::sd(successful$specificity, na.rm = TRUE),
      mean_precision = mean(successful$precision, na.rm = TRUE),
      sd_precision = stats::sd(successful$precision, na.rm = TRUE),
      mean_fdr = mean(successful$false_discovery_rate, na.rm = TRUE),
      sd_fdr = stats::sd(successful$false_discovery_rate, na.rm = TRUE),
      mean_projection_distance = mean(successful$projection_distance, na.rm = TRUE),
      sd_projection_distance = stats::sd(successful$projection_distance, na.rm = TRUE),
      mean_c_index = mean(successful$c_index, na.rm = TRUE),
      sd_c_index = stats::sd(successful$c_index, na.rm = TRUE),
      mean_runtime_seconds = mean(successful$runtime_seconds, na.rm = TRUE),
      sd_runtime_seconds = stats::sd(successful$runtime_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    method_selection_table <- as.data.frame(
      table(
        covariance = successful$best_covariance,
        d = successful$best_d
      )
    )

    names(method_selection_table) <- c(
      "covariance",
      "d",
      "frequency"
    )
  } else {
    summary_table <- data.frame()
    method_selection_table <- data.frame()
  }

  overall_end <- Sys.time()

  total_runtime_seconds <- as.numeric(
    difftime(
      overall_end,
      overall_start,
      units = "secs"
    )
  )

  out <- list(
    replication_results = replication_results,
    results = result_table,
    summary = summary_table,
    method_selection = method_selection_table,
    settings = list(
      n = n,
      p = p,
      active = active,
      dimension_regime = dimension_regime,
      covariance = covariance,
      rho = rho,
      signal_strength = signal_strength,
      censoring_rate = censoring_rate,
      n_rep = n_rep,
      screening = screening,
      screen_top = screen_top,
      covariance_methods = covariance_methods,
      d_grid = d_grid,
      n_slices = n_slices,
      selector = selector,
      variable_top = variable_top,
      seed = seed
    ),
    total_runtime_seconds = total_runtime_seconds
  )

  class(out) <- "ssurvdr_simulation_run"

  out
}
