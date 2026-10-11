# ssurvdr

Survival-informed SAVE dimension reduction with sample or maximum-entropy covariance stabilisation, reduced Cox modelling and information-criterion selection.

This is release candidate 0.1.0, prepared for initial CRAN submission.

## Small reproducible example

```r
library(ssurvdr)
dat <- simulate_ssurvdr(n = 100, p = 12, active = 1:3, seed = 2026)
fit <- fit_ssurvdr(dat$X, dat$time, dat$status, d = 2,
                  n_slices = 3, covariance = "mec", backend = "primal")
fit$criteria
path <- fit_ssurvdr_sparse_path(dat$X, dat$time, dat$status,
                               dense_fit = fit, support_grid = c(4, 8, 12))
path$path
```

## Statistical interpretation

- `fit_ssurvdr()` estimates a dense basis; its name does not imply that it automatically selects a sparse support.
- `fit_ssurvdr_sparse_path()` ranks rows of a fitted basis and refits the reduced Cox model at candidate support sizes, including full support. This is post-estimation sparsification.
- `sparse_projection()` thresholds a dense primal SAVE kernel and does not support the compact dual kernel.
- `compute_icomp()` calculates C1 covariance complexity from stabilised eigenvalues and uses it in reduced Cox information criteria. No original-predictor support-size penalty is added.
- `compute_c1f_weights()` retains its historical name and `c1f` field. That diagnostic uses the implemented C1 arithmetic/geometric eigenvalue formula. The weights depend on estimated basis row norms and a separate eigenvalue-dispersion adjustment.
- `evaluate_ssurvdr()` reports the fitted training concordance. Independent predictive assessment requires a held-out sample or cross-validation with all preprocessing, screening and selection performed inside training folds.

The research scripts accompany the audit package outside the installable source. They are not run by package checks.
