#' Compute information criteria for reduced Cox model
#'
#' Computes AIC, BIC, CAIC, ICOMP(IFIM), and ICOMP(misspec) for a reduced
#' Cox proportional hazards model fitted after sufficient dimension reduction.
#'
#' @param cox_fit Object returned by [fit_reduced_cox()].
#' @param eps Small positive stabilisation constant for eigenvalues.
#'
#' @return A list containing information criteria, covariance-complexity
#'   diagnostics, and likelihood quantities.
#'
#' @export
compute_icomp <- function(
    cox_fit,
    eps = 1e-8
) {
  if (missing(cox_fit)) {
    stop("'cox_fit' must be supplied.", call. = FALSE)
  }

  if (!inherits(cox_fit, "ssurvdr_cox")) {
    stop(
      "'cox_fit' must be an object returned by fit_reduced_cox().",
      call. = FALSE
    )
  }

  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("'eps' must be a single positive numeric value.", call. = FALSE)
  }

  model <- cox_fit$cox_model

  loglik <- as.numeric(cox_fit$loglik_final)

  if (!is.finite(loglik)) {
    stop("The fitted Cox model does not contain a finite log-likelihood.", call. = FALSE)
  }

  n <- cox_fit$n
  k <- length(stats::coef(model))

  if (k < 1L) {
    stop("The fitted Cox model must contain at least one coefficient.", call. = FALSE)
  }

  covariance_complexity <- function(V) {
    V <- as.matrix(V)
    V <- (V + t(V)) / 2

    eigen_values <- eigen(
      V,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    eigen_values <- pmax(eigen_values, eps)

    s <- length(eigen_values)

    arithmetic_mean <- mean(eigen_values)

    geometric_mean <- exp(mean(log(eigen_values)))

    c1 <- (s / 2) * log(arithmetic_mean / geometric_mean)

    list(
      c1 = as.numeric(c1),
      eigenvalues = eigen_values,
      arithmetic_mean = arithmetic_mean,
      geometric_mean = geometric_mean,
      condition_number = max(eigen_values) / min(eigen_values)
    )
  }

  vcov_ifim <- stats::vcov(model)

  if (anyNA(vcov_ifim) || any(!is.finite(vcov_ifim))) {
    stop("The IFIM covariance matrix contains non-finite values.", call. = FALSE)
  }

  ifim_complexity <- covariance_complexity(vcov_ifim)

  score_residuals <- tryCatch(
    {
      stats::residuals(
        model,
        type = "score"
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(score_residuals)) {
    vcov_misspec <- matrix(
      NA_real_,
      nrow = k,
      ncol = k
    )

    misspec_complexity <- list(
      c1 = NA_real_,
      eigenvalues = rep(NA_real_, k),
      arithmetic_mean = NA_real_,
      geometric_mean = NA_real_,
      condition_number = NA_real_
    )
  } else {
    score_residuals <- as.matrix(score_residuals)

    if (ncol(score_residuals) != k) {
      score_residuals <- matrix(
        score_residuals,
        ncol = k
      )
    }

    G <- t(score_residuals) %*% score_residuals

    vcov_misspec <- vcov_ifim %*% G %*% vcov_ifim
    vcov_misspec <- (vcov_misspec + t(vcov_misspec)) / 2

    misspec_complexity <- covariance_complexity(vcov_misspec)
  }

  aic <- -2 * loglik + 2 * k

  bic <- -2 * loglik + k * log(n)

  caic <- -2 * loglik + k * (log(n) + 1)

  icomp_ifim <- -2 * loglik + 2 * ifim_complexity$c1

  icomp_misspec <- -2 * loglik + 2 * misspec_complexity$c1

  criteria <- data.frame(
    criterion = c(
      "AIC",
      "BIC",
      "CAIC",
      "ICOMP_IFIM",
      "ICOMP_MISSPEC"
    ),
    value = c(
      aic,
      bic,
      caic,
      icomp_ifim,
      icomp_misspec
    ),
    stringsAsFactors = FALSE
  )

  out <- list(
    criteria = criteria,
    aic = aic,
    bic = bic,
    caic = caic,
    icomp_ifim = icomp_ifim,
    icomp_misspec = icomp_misspec,
    loglik = loglik,
    n = n,
    k = k,
    vcov_ifim = vcov_ifim,
    vcov_misspec = vcov_misspec,
    c1_ifim = ifim_complexity$c1,
    c1_misspec = misspec_complexity$c1,
    ifim_diagnostics = ifim_complexity,
    misspec_diagnostics = misspec_complexity,
    cox_fit = cox_fit
  )

  class(out) <- "ssurvdr_icomp"

  out
}
