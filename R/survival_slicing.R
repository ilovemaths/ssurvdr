#' Construct survival-informed slices for censored survival data
#'
#' This function partitions survival observations by censoring status and then
#' slices observed survival times within each censoring group. It is intended
#' for sufficient dimension reduction methods for censored survival data.
#'
#' @param time A numeric vector of observed survival times.
#' @param status A numeric, integer, or logical vector of censoring indicators.
#'   Conventionally, 1 indicates observed event and 0 indicates censoring.
#' @param n_slices Number of time slices to construct within each censoring group.
#' @param min_slice_size Minimum number of observations allowed in a slice.
#' @param drop_empty Logical. If TRUE, empty slices are dropped.
#'
#' @return A list with the slice membership, slice labels, slice table, and
#'   basic metadata.
#'
#' @export
survival_slicing <- function(
    time,
    status,
    n_slices = 5,
    min_slice_size = 2,
    drop_empty = TRUE
) {
  if (missing(time) || missing(status)) {
    stop("Both 'time' and 'status' must be supplied.", call. = FALSE)
  }

  if (!is.numeric(time)) {
    stop("'time' must be a numeric vector.", call. = FALSE)
  }

  if (anyNA(time)) {
    stop("'time' must not contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(time))) {
    stop("'time' must contain only finite values.", call. = FALSE)
  }

  if (any(time < 0)) {
    stop("'time' must not contain negative survival times.", call. = FALSE)
  }

  if (length(time) != length(status)) {
    stop("'time' and 'status' must have the same length.", call. = FALSE)
  }

  if (length(time) < 2) {
    stop("At least two observations are required.", call. = FALSE)
  }

  if (anyNA(status)) {
    stop("'status' must not contain missing values.", call. = FALSE)
  }

  if (is.logical(status)) {
    status <- as.integer(status)
  }

  if (!is.numeric(status) && !is.integer(status)) {
    stop("'status' must be numeric, integer, or logical.", call. = FALSE)
  }

  status_values <- sort(unique(as.integer(status)))

  if (!all(status_values %in% c(0L, 1L))) {
    stop("'status' must contain only 0/1 values, where 1 = event and 0 = censored.", call. = FALSE)
  }

  status <- as.integer(status)

  if (!is.numeric(n_slices) || length(n_slices) != 1L || n_slices < 1) {
    stop("'n_slices' must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(min_slice_size) || length(min_slice_size) != 1L || min_slice_size < 1) {
    stop("'min_slice_size' must be a positive integer.", call. = FALSE)
  }

  n_slices <- as.integer(n_slices)
  min_slice_size <- as.integer(min_slice_size)

  n <- length(time)

  slice_id <- rep(NA_integer_, n)
  slice_label <- rep(NA_character_, n)

  slice_counter <- 0L
  slice_records <- list()

  for (current_status in c(0L, 1L)) {
    idx <- which(status == current_status)

    if (length(idx) == 0L) {
      next
    }

    group_n <- length(idx)

    group_slices <- min(n_slices, ceiling(group_n / min_slice_size))

    if (group_slices < 1L) {
      group_slices <- 1L
    }

    ordered_idx <- idx[order(time[idx], seq_along(idx))]

    cut_groups <- cut(
      seq_along(ordered_idx),
      breaks = group_slices,
      labels = FALSE,
      include.lowest = TRUE
    )

    for (h in seq_len(group_slices)) {
      members <- ordered_idx[cut_groups == h]

      if (length(members) == 0L && isTRUE(drop_empty)) {
        next
      }

      slice_counter <- slice_counter + 1L

      current_label <- paste0(
        ifelse(current_status == 1L, "event", "censored"),
        "_slice_",
        h
      )

      slice_id[members] <- slice_counter
      slice_label[members] <- current_label

      slice_records[[length(slice_records) + 1L]] <- data.frame(
        slice_id = slice_counter,
        status = current_status,
        slice_within_status = h,
        label = current_label,
        n = length(members),
        time_min = if (length(members) > 0L) min(time[members]) else NA_real_,
        time_max = if (length(members) > 0L) max(time[members]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (anyNA(slice_id)) {
    stop("Internal slicing error: some observations were not assigned to slices.", call. = FALSE)
  }

  slice_table <- do.call(rbind, slice_records)
  rownames(slice_table) <- NULL

  out <- list(
    slice_id = slice_id,
    slice_label = slice_label,
    slice_table = slice_table,
    n = n,
    n_slices_requested = n_slices,
    n_slices_created = length(unique(slice_id)),
    min_slice_size = min_slice_size,
    time = time,
    status = status
  )

  class(out) <- "ssurvdr_slices"

  out
}
