## Manual checks for unxar1() and homcsxar1().
##
## Run from the repository root with:
##   Rscript notes/check_xar1_covstruct_v001.R
##
## The script loads the local glmmTMB source tree, exercises the new
## covariance structures, and prints the output most useful for inspection.

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

section <- function(x) {
    cat("\n", paste(rep("=", nchar(x) + 4), collapse = ""), "\n", sep = "")
    cat("  ", x, "\n", sep = "")
    cat(paste(rep("=", nchar(x) + 4), collapse = ""), "\n\n", sep = "")
}

rho_to_theta <- function(rho) rho / sqrt(1 - rho^2)
homcs_to_theta <- function(rho, n) {
    a <- 1 / (n - 1)
    qlogis((rho + a) / (1 + a))
}

summarize_xar1 <- function(fit) {
    vc <- VarCorr(fit)$cond[[1]]
    cc <- attr(vc, "correlation")
    ss <- attr(vc, "stddev")
    data.frame(
        sd_member_1 = ss[[1]],
        sd_member_2 = ss[[2]],
        rho_member = cc[1, 2],
        phi = cc[1, 3],
        cross_lag1 = cc[1, 4],
        row.names = NULL
    )
}

make_dyad <- function(ng = 80, nt = 5, rho = 0.35, phi = 0.55,
                      sd = c(0.8, 1.2), beta = 1) {
    member <- 1:2
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_time <- outer(seq_len(nt), seq_len(nt), function(i, j) phi ^ abs(i - j))
    R <- kronecker(R_time, R_member)
    sd_coord <- rep(sd, nt)
    Sigma <- R * outer(sd_coord, sd_coord)
    L <- chol(Sigma)

    dd <- expand.grid(member = member, time = seq_len(nt), group = factor(seq_len(ng)))
    u <- do.call(c, replicate(ng, beta + as.numeric(t(L) %*% rnorm(2 * nt)), simplify = FALSE))
    dd$y <- u
    dd
}

section("1. Validation errors for malformed coordinate factors")
bad <- expand.grid(group = factor(1:3), member = 1:2, time = 1:3)
bad$y <- rnorm(nrow(bad))
cat("Intercept error:\n")
print(try(glmmTMB(y ~ 1 + homcsxar1(membertime(member, time) | group),
                  data = bad, doFit = FALSE), silent = TRUE))

bad$time_only <- numFactor(bad$time)
cat("\nOne-dimensional numFactor error:\n")
print(try(glmmTMB(y ~ 1 + homcsxar1(time_only + 0 | group),
                  data = bad, doFit = FALSE), silent = TRUE))

section("2. Simulated Gaussian dyad: unxar1")
set.seed(11)
simdat <- make_dyad()
fit_un <- glmmTMB(
    y ~ 1 + unxar1(membertime(member=member, time) + 0 | group),
    data = simdat,
    dispformula = ~0
)
print(fit_un)
print(summary(fit_un))
print(VarCorr(fit_un), maxdim = 10)
cat("\nCompact covariance summary:\n")
print(summarize_xar1(fit_un))
cat("\nWald confidence intervals:\n")
print(confint(fit_un, method = "wald"))

section("3. Missing observations within groups")
missdat <- simdat[!(simdat$group %in% c("1", "2") &
                    simdat$member == 2 & simdat$time %in% c(4, 5)), ]
full_grid <- expand.grid(member = 1:2, time = sort(unique(simdat$time)))
full_member_time <- numFactor(full_grid$member, full_grid$time)
obs_member_time <- numFactor(missdat$member, missdat$time)
missdat$member_time <- factor(as.character(obs_member_time),
                              levels = levels(full_member_time))

fit_miss <- glmmTMB(
    y ~ 1 + homcsxar1(member_time + 0 | group),
    data = missdat,
    dispformula = ~0
)
print(fit_miss)
print(VarCorr(fit_miss), maxdim = 10)
cat("\nCompact covariance summary:\n")
print(summarize_xar1(fit_miss))

section("4. Real-data stress check: Salamanders")
data(Salamanders)
Salamanders$spp_num <- as.numeric(Salamanders$spp)
Salamanders$sample_num <- as.numeric(Salamanders$sample)
Salamanders$spp_sample <- numFactor(Salamanders$spp_num, Salamanders$sample_num)

cat("Grid dimensions:\n")
cat("sites =", length(unique(Salamanders$site)),
    "species =", length(unique(Salamanders$spp)),
    "samples =", length(unique(Salamanders$sample)),
    "rows =", nrow(Salamanders), "\n")

fit_sal <- glmmTMB(
    count ~ mined + homcsxar1(spp_sample + 0 | site),
    data = Salamanders,
    family = poisson,
    se = FALSE,
    control = glmmTMBControl(optCtrl = list(iter.max = 100, eval.max = 150))
)
print(fit_sal)
print(VarCorr(fit_sal), maxdim = 8)
cat("\nFirst 8x8 block of the reported correlation matrix:\n")
print(round(attr(VarCorr(fit_sal)$cond[[1]], "correlation")[1:8, 1:8], 3))

section("Done")
cat("If these fits converge and the compact summaries look plausible, the next\n")
cat("step is a replicated simulation-recovery script for bias/RMSE/coverage.\n")
