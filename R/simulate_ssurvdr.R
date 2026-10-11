#' Simulate high-dimensional censored survival data for ssurvdr
#'
#' Simulates predictors, sparse true SDR direction, event times, censoring times,
#' observed survival times, and censoring indicators for testing sparse survival
#' dimension reduction methods.
#'
#' @param n Sample size.
#' @param p Number of predictors.
#' @param active Integer vector of active predictor indices.
#' @param beta Optional numeric vector of active coefficients. If NULL, equal
#'   coefficients are used.
#' @param covariance Covariance structure. One of "identity", "ar1",
#'   "compound", or "block".
#' @param rho Correlation parameter for AR(1), compound, or block covariance.
#' @param block_size Block size for block covariance.
#' @param signal_strength Multiplicative signal strength.
#' @param baseline_hazard Exponential baseline hazard rate.
#' @param censoring_rate Rate parameter for exponential censoring distribution.
#' @param standardise_X Logical. If TRUE, predictors are centred and scaled.
#' @param seed Optional random seed.
#'
#' @return A list containing simulated predictors, survival outcome, true basis,
#'   true active set, covariance matrix, and simulation metadata.
#'
#' @export
simulate_ssurvdr <- function(
    n = 100,
    p = 50,
    active = 1:5,
    beta = NULL,
    covariance = c("identity", "ar1", "compound", "block"),
    rho = 0.5,
    block_size = 10,
    signal_strength = 1,
    baseline_hazard = 0.1,
    censoring_rate = 0.1,
    standardise_X = TRUE,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  covariance <- match.arg(covariance)

  if (!is.numeric(n) || length(n) != 1L || n < 2) {
    stop("'n' must be a single integer greater than or equal to 2.", call. = FALSE)
  }

  if (!is.numeric(p) || length(p) != 1L || p < 1) {
    stop("'p' must be a single positive integer.", call. = FALSE)
  }

  n <- as.integer(n)
  p <- as.integer(p)

  if (!is.numeric(active) || length(active) < 1L) {
    stop("'active' must be a non-empty integer vector.", call. = FALSE)
  }

  active <- sort(unique(as.integer(active)))

  if (any(active < 1L) || any(active > p)) {
    stop("All entries of 'active' must be valid predictor indices.", call. = FALSE)
  }

  s <- length(active)

  if (is.null(beta)) {
    beta <- rep(1 / sqrt(s), s)
  }

  if (!is.numeric(beta) || length(beta) != s) {
    stop("'beta' must be NULL or have the same length as 'active'.", call. = FALSE)
  }

  if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1) {
    stop("'rho' must be a single numeric value with absolute value less than 1.", call. = FALSE)
  }

  if (!is.numeric(block_size) || length(block_size) != 1L || block_size < 1) {
    stop("'block_size' must be a positive integer.", call. = FALSE)
  }

  block_size <- as.integer(block_size)

  if (!is.numeric(signal_strength) || length(signal_strength) != 1L) {
    stop("'signal_strength' must be a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(baseline_hazard) ||
      length(baseline_hazard) != 1L ||
      baseline_hazard <= 0) {
    stop("'baseline_hazard' must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(censoring_rate) ||
      length(censoring_rate) != 1L ||
      censoring_rate <= 0) {
    stop("'censoring_rate' must be a single positive numeric value.", call. = FALSE)
  }

  Sigma <- diag(p)

  if (covariance == "ar1") {
    idx <- seq_len(p)

    Sigma <- outer(
      idx,
      idx,
      function(j, k) rho^abs(j - k)
    )
  }

  if (covariance == "compound") {
    Sigma <- matrix(rho, nrow = p, ncol = p)
    diag(Sigma) <- 1
  }

  if (covariance == "block") {
    Sigma <- diag(p)

    starts <- seq(1L, p, by = block_size)

    for (start in starts) {
      end <- min(start + block_size - 1L, p)
      block_index <- start:end

      Sigma[block_index, block_index] <- rho
      diag(Sigma)[block_index] <- 1
    }
  }

  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma
  )

  X <- as.matrix(X)

  if (isTRUE(standardise_X)) {
    X <- scale(X, center = TRUE, scale = TRUE)
    X <- as.matrix(X)
  }

  true_basis <- rep(0, p)
  true_basis[active] <- beta
  true_basis <- true_basis / sqrt(sum(true_basis^2))

  eta <- as.numeric(signal_strength * X %*% true_basis)

  U <- stats::runif(n)

  true_time <- -log(U) / (baseline_hazard * exp(eta))

  censoring_time <- stats::rexp(
    n = n,
    rate = censoring_rate
  )

  time <- pmin(true_time, censoring_time)

  status <- as.integer(true_time <= censoring_time)

  censoring_proportion <- mean(status == 0L)

  colnames(X) <- paste0("X", seq_len(p))

  out <- list(
    X = X,
    time = time,
    status = status,
    true_time = true_time,
    censoring_time = censoring_time,
    true_basis = matrix(true_basis, ncol = 1),
    active = active,
    beta = beta,
    eta = eta,
    Sigma = Sigma,
    censoring_proportion = censoring_proportion,
    metadata = list(
      n = n,
      p = p,
      active = active,
      covariance = covariance,
      rho = rho,
      block_size = block_size,
      signal_strength = signal_strength,
      baseline_hazard = baseline_hazard,
      censoring_rate = censoring_rate,
      standardise_X = standardise_X
    )
  )

  class(out) <- "ssurvdr_simulation"

  out
}
