## Adaptation of notes/kronecker.qmd for the new unxar1() structure.
##
## kronecker.qmd simulates AR1(row) x AR1(col).  v0.0.1 does not implement a
## parsimonious AR1 x AR1 covariance, but unxar1() can represent the same
## covariance as UN(row) x AR1(col).  This script checks that recovery.
##
## Run from the repository root:
##   Rscript notes/check_kronecker_qmd_unxar1.R
##
## Optional full stress-test dimensions from kronecker.qmd:
##   FULL_KRON_CHECK=1 Rscript notes/check_kronecker_qmd_unxar1.R

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo_root <- if (length(cmd_file)) {
    normalizePath(file.path(dirname(sub("^--file=", "", cmd_file[[1]])), ".."))
} else {
    normalizePath(getwd())
}
pkg_dir <- file.path(repo_root, "glmmTMB")

if (requireNamespace("pkgload", quietly = TRUE) && dir.exists(pkg_dir)) {
    pkgload::load_all(pkg_dir, quiet = TRUE, helpers = FALSE)
} else {
    library(glmmTMB)
}

rho_to_theta <- function(rho) rho / sqrt(1 - rho^2)

fit_kron_check <- function(n_R = 5, n_C = 6, n_block = 40, n_rep = 3,
                           rho_R = 0.8, rho_C = 0.2,
                           sigma_RE = 2, sigma_eps = 1,
                           seed = 101) {
    set.seed(seed)
    R_row <- outer(seq_len(n_R), seq_len(n_R), function(i, j) rho_R^abs(i - j))
    R_col <- outer(seq_len(n_C), seq_len(n_C), function(i, j) rho_C^abs(i - j))
    Sigma <- kronecker(R_col, sigma_RE^2 * R_row)
    L <- chol(Sigma)
    U <- t(replicate(n_block, as.numeric(t(L) %*% rnorm(n_R * n_C))))

    dd <- expand.grid(row = seq_len(n_R), col = seq_len(n_C),
                      block = factor(seq_len(n_block)), rep = seq_len(n_rep))
    coord <- dd$row + (dd$col - 1) * n_R
    dd$u <- U[cbind(as.integer(dd$block), coord)]
    dd$y <- dd$u + rnorm(nrow(dd), sd = sigma_eps)

    start_theta <- c(rep(log(sigma_RE), n_R), put_cor(R_row), rho_to_theta(rho_C))
    fit <- glmmTMB(
        y ~ 1 + unxar1(membertime(row, col) + 0 | block),
        data = dd,
        REML = TRUE,
        start = list(beta = 0, betadisp = log(sigma_eps), theta = start_theta),
        control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1200))
    )

    vc <- VarCorr(fit)$cond[[1]]
    cc <- attr(vc, "correlation")
    ss <- attr(vc, "stddev")
    row_lag1 <- cc[cbind(1:(n_R - 1), 2:n_R)]
    row_lag2 <- cc[cbind(1:(n_R - 2), 3:n_R)]

    list(
        fit = fit,
        summary = data.frame(
            n_R = n_R,
            n_C = n_C,
            n_block = n_block,
            n_rep = n_rep,
            convergence = fit$fit$convergence,
            pdHess = fit$sdr$pdHess,
            true_sigma_RE = sigma_RE,
            est_sd_mean = mean(ss),
            est_sd_range_low = min(ss),
            est_sd_range_high = max(ss),
            true_rho_row_lag1 = rho_R,
            est_rho_row_lag1_mean = mean(row_lag1),
            true_rho_row_lag2 = rho_R^2,
            est_rho_row_lag2_mean = mean(row_lag2),
            true_rho_col_lag1 = rho_C,
            est_rho_col_lag1 = cc[1, 1 + n_R],
            true_cross_row_col = rho_R * rho_C,
            est_cross_row_col = cc[1, 2 + n_R],
            true_sigma_eps = sigma_eps,
            est_sigma_eps = sigma(fit),
            theta_length = length(getME(fit, "theta")),
            row.names = NULL
        )
    )
}

cat("\nScaled-down kronecker.qmd check: AR1(row) x AR1(col) fitted as UN(row) x AR1(col)\n\n")
res <- fit_kron_check()
print(res$summary)

if (identical(Sys.getenv("FULL_KRON_CHECK"), "1")) {
    cat("\nFull kronecker.qmd dimensions: 10 x 10 x 10 blocks x 3 reps\n")
    cat("This is a dense 100-dimensional random-effect block with 56 theta parameters and may take a while.\n\n")
    full <- fit_kron_check(n_R = 10, n_C = 10, n_block = 10, n_rep = 3)
    print(full$summary)
}
