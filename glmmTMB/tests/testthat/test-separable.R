stopifnot(require("testthat"), require("glmmTMB"))

make_sep_dat <- function(n_member = 2, n_time = 3, n_group = 2,
                         reps = FALSE) {
    args <- list(member = factor(paste0("m", seq_len(n_member))),
                 time = factor(seq_len(n_time)),
                 group = factor(seq_len(n_group)))
    if (reps) args$rep <- 1:2
    dd <- do.call(expand.grid, args)
    dd$y <- seq_len(nrow(dd)) / nrow(dd)
    dd
}

ar1_to_theta <- function(phi) phi / sqrt(1 - phi^2)

homcs_to_theta <- function(rho, n) {
    lower <- -1 / (n - 1)
    qlogis((rho - lower) / (1 - lower))
}

fit_fixed_theta <- function(form, dd, theta) {
    glmmTMB(form, data = dd,
            start = list(theta = theta),
            map = list(theta = factor(rep(NA, length(theta)))))
}

joint_nll_at <- function(form, dd, theta, b, sigma = 2) {
    fit <- glmmTMB(form, data = dd,
                   start = list(beta = 0, betadisp = log(sigma),
                                theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))
    p <- fit$obj$env$last.par.best
    p[names(p) == "beta"] <- 0
    p[names(p) == "betadisp"] <- log(sigma)
    p[names(p) == "b"] <- b
    p[names(p) == "theta"] <- theta
    fit$obj$env$f(p)
}

ar1_corr <- function(n, phi) {
    outer(seq_len(n), seq_len(n), function(i, j) phi^abs(i - j))
}

cs_corr <- function(n, rho) {
    R <- matrix(rho, n, n)
    diag(R) <- 1
    R
}

toep_corr <- function(rho) {
    n <- length(rho) + 1L
    R <- diag(n)
    for (lag in seq_len(n - 1L)) {
        R[row(R) == col(R) + lag] <- rho[[lag]]
        R[col(R) == row(R) + lag] <- rho[[lag]]
    }
    R
}

check_dense_reference <- function(form, dense_form, dd, theta, sd_full,
                                  R_full, check_vc = TRUE) {
    dd$y <- 0
    theta_dense <- c(log(sd_full), put_cor(R_full))
    b <- seq(-0.4, 0.5, length.out = length(sd_full) * nlevels(dd$group))

    sep_nll <- joint_nll_at(form, dd, theta, b)
    dense_nll <- joint_nll_at(dense_form, dd, theta_dense, b)
    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)

    if (check_vc) {
        fit <- fit_fixed_theta(form, dd, theta)
        vc <- VarCorr(fit)$cond[[1]]
        expect_equal(unname(attr(vc, "stddev")), sd_full, tolerance = 1e-6)
        expect_equal(unname(attr(vc, "correlation")), R_full, tolerance = 1e-6)
    }
}

test_that("sepgrid builds complete two-dimensional levels", {
    member <- factor(c("A", "B"), levels = c("A", "B"))
    time <- factor(c(1, 3), levels = 1:3)
    grid <- sepgrid(member, time)

    expect_equal(nlevels(grid), 6)
    expect_equal(unname(parseNumLevels(levels(grid))),
                 unname(as.matrix(expand.grid(1:2, 1:3))),
                 check.attributes = FALSE)
})

test_that("glmmTMB preserves separable levels without preserving unrelated levels", {
    dd <- expand.grid(member = factor(c("A", "B")),
                      time = factor(c(1, 3), levels = 1:3),
                      group = factor(1:2),
                      rep = 1:2)
    dd$fixed_factor <- factor(rep(c("a", "b"), length.out = nrow(dd)),
                              levels = c("a", "b", "unused"))
    dd$y <- 0

    fit <- glmmTMB(y ~ fixed_factor +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepDims, c(2L, 3L))
    expect_equal(levels(fit$fr$time), as.character(1:3))
    expect_equal(levels(fit$fr$fixed_factor), c("a", "b"))
})

test_that("separable parser flattens product chains and records scale syntax", {
    f <- y ~ 1 +
        separable(us(0 + member) %x% ar1(0 + time) %x% cs(0 + item) | group,
                  scale = product(us(0 + member), cs(0 + item)))
    ss <- reformulas::splitForm(f, specials = c(names(.valid_covstruct), "s"))
    spec <- glmmTMB:::.sep_specs_from_split(ss)[[1]]

    expect_equal(spec$grid, c("member", "time", "item"))
    expect_equal(spec$margins$struc, c("us", "ar1", "cs"))
    expect_equal(spec$scale$mode, "selected_product")
    expect_equal(spec$scale$margins$struc, c("us", "cs"))
})

test_that("separable dense x ar1 matches dense MVN", {
    dd <- make_sep_dat()
    rho <- 0.3
    phi <- 0.45
    sd_member <- c(0.8, 1.2)
    R_member <- cs_corr(2, rho)
    R_time <- ar1_corr(3, phi)
    sd_full <- rep(sd_member, 3)
    R_full <- kronecker(R_time, R_member)

    check_dense_reference(
        y ~ 1 + separable(cs(0 + member) %x% ar1(0 + time) | group),
        y ~ 1 + us(sepgrid(member, time) + 0 | group),
        dd,
        c(log(sd_member), homcs_to_theta(rho, 2), ar1_to_theta(phi)),
        sd_full, R_full
    )
})

test_that("separable product scale matches dense MVN", {
    dd <- expand.grid(member = factor(paste0("m", 1:2)),
                      item = factor(paste0("i", 1:2)),
                      group = factor(1:2))
    dd$y <- 0
    sd_member <- c(0.7, 1.1)
    sd_item <- c(1.2, 1.4)
    R_member <- diag(2)
    R_item <- cs_corr(2, 0.25)
    sd_full <- as.vector(outer(sd_member, sd_item))
    R_full <- kronecker(R_item, R_member)

    check_dense_reference(
        y ~ 1 + separable(diag(0 + member) %x% us(0 + item) | group,
                          scale = product()),
        y ~ 1 + us(sepgrid(member, item) + 0 | group),
        dd,
        c(log(sd_member), log(sd_item), put_cor(R_item)),
        sd_full, R_full
    )
})

test_that("separable supports Toeplitz margins", {
    dd <- expand.grid(time = factor(1:3),
                      member = factor(paste0("m", 1:2)),
                      group = factor(1:2))
    dd$y <- 0
    rho <- c(0.25, 0.1)
    phi <- 0.35
    sd_time <- c(0.8, 1.0, 1.2)
    R_time <- toep_corr(rho)
    R_member <- ar1_corr(2, phi)
    sd_full <- rep(sd_time, 2)
    R_full <- kronecker(R_member, R_time)

    check_dense_reference(
        y ~ 1 + separable(toep(0 + time) %x% ar1(0 + member) | group),
        y ~ 1 + us(sepgrid(time, member) + 0 | group),
        dd,
        c(log(sd_time), ar1_to_theta(rho), ar1_to_theta(phi)),
        sd_full, R_full
    )
})

test_that("separable supports fixed covariance and spatial margins", {
    dd <- expand.grid(member = factor(c("a", "b")),
                      time = numFactor(c(0, 1, 3)),
                      group = factor(1:2))
    dd$y <- 0
    K <- matrix(c(1, 0.2, 0.2, 1.4), 2, 2,
                dimnames = list(levels(dd$member), levels(dd$member)))
    global_sd <- 1.25
    decay <- -0.2
    R_member <- cov2cor(K)
    R_time <- exp(-exp(decay) * as.matrix(stats::dist(c(0, 1, 3))))
    diag(R_time) <- 1
    sd_full <- unname(global_sd * rep(sqrt(diag(K)), 3))
    R_full <- kronecker(R_time, R_member)
    env <- list2env(list(K = K), parent = environment())
    form <- y ~ 1 +
        separable(propto(0 + member, K) %x% ou(0 + time) | group,
                  scale = global())
    dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
    environment(form) <- env
    environment(dense_form) <- env

    fit <- glmmTMB(form, data = dd, doFit = FALSE)
    restruc <- fit$condReStruc[[1]]
    expect_equal(restruc$sepMatrixPayloadKinds, c(1L, 2L))

    check_dense_reference(form, dense_form, dd, c(log(global_sd), decay),
                          sd_full, R_full)
})

test_that("separable supports fixed covariance without estimated scale", {
    dd <- expand.grid(member = factor(c("a", "b")),
                      time = factor(1:3),
                      group = factor(1:2))
    dd$y <- 0
    K <- matrix(c(1, 0.2, 0.2, 1.4), 2, 2)
    phi <- 0.35
    R_member <- cov2cor(K)
    R_time <- ar1_corr(3, phi)
    sd_full <- unname(rep(sqrt(diag(K)), 3))
    R_full <- kronecker(R_time, R_member)
    env <- list2env(list(K = K), parent = environment())
    form <- y ~ 1 +
        separable(equalto(0 + member, K) %x% ar1(0 + time) | group)
    dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
    environment(form) <- env
    environment(dense_form) <- env

    check_dense_reference(form, dense_form, dd, ar1_to_theta(phi),
                          sd_full, R_full)
})

test_that("separable supports three-margin products", {
    dd <- expand.grid(member = factor(c("a", "b")),
                      time = factor(1:2),
                      item = factor(c("i1", "i2")),
                      group = factor(1:2))
    dd$y <- 0

    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) %x%
                                     homtoep(0 + item) | group,
                                 scale = global()),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepDims, c(2L, 2L, 2L))
    expect_equal(fit$condReStruc[[1]]$blockNumTheta, 4L)
})

test_that("separable product syntax supports multi-column dense margins", {
    dd <- make_sep_dat()
    dd$x <- seq_len(nrow(dd))

    fit <- glmmTMB(y ~ 1 +
                       separable(us(0 + member + member:x) %x%
                                     ar1(0 + time) | group,
                                 scale = us(0 + member + member:x)),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepDims, c(4L, 3L))
    expect_equal(unname(fit$condReStruc[[1]]$blockSize), 12)
    expect_s3_class(fit$condList$reXterms[[1]], "separable_reXterms")
})

test_that("multiple separable terms and smooth augmentation keep spec order", {
    skip_if_not_installed("mgcv")
    s <- mgcv::s
    dd <- make_sep_dat()
    dd$x <- seq_len(nrow(dd)) / nrow(dd)

    fit <- glmmTMB(y ~ s(x, k = 4) +
                       separable(us(0 + member) %x% ar1(0 + time) | group) +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condList$ss, c("homdiag", "separable", "separable"))
    expect_equal(fit$condReStruc[[2]]$sepDims, c(2L, 3L))
    expect_equal(fit$condReStruc[[3]]$sepDims, c(2L, 3L))
})

test_that("separable preserves user-facing formulas", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1 + separable(homcs(0 + member) %x%
                                         ar1(0 + time) | group),
                   ziformula = ~ separable(homcs(0 + member) %x%
                                               ar1(0 + time) | group),
                   dispformula = ~ separable(homcs(0 + member) %x%
                                                 ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    txt <- paste(deparse(fit$call$formula), collapse = " ")
    ztxt <- paste(deparse(fit$call$ziformula), collapse = " ")
    dtxt <- paste(deparse(fit$call$dispformula), collapse = " ")
    expect_match(txt, "homcs\\(0 \\+ member\\)")
    expect_match(ztxt, "homcs\\(0 \\+ member\\)")
    expect_match(dtxt, "homcs\\(0 \\+ member\\)")
    expect_false(grepl("data.frame", paste(txt, ztxt, dtxt), fixed = TRUE))
})

test_that("separable validates syntax and scale choices", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 + separable(us(member) %x% ar1(0 + time) | group,
                                  scale = us(member)),
                data = dd, doFit = FALSE),
        "no-intercept"
    )
    expect_error(
        glmmTMB(y ~ 1 + separable(us(0 + member) %x% homcs(0 + time) | group),
                data = dd, doFit = FALSE),
        "specify the scale mode"
    )
    expect_error(
        glmmTMB(y ~ 1 + separable(ar1(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "no unambiguous scale margin"
    )
    expect_error(
        glmmTMB(y ~ 1 + separable(foo(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) margin: foo"
    )
    expect_error(
        glmmTMB(y ~ 1 + separable(us(0 + member) %x% ar1(0 + time) | group,
                                  scale = product(cs(0 + item))),
                data = dd, doFit = FALSE),
        "must match one of the specified margins"
    )
})

test_that("separable dense x ar1 models fit successfully", {
    set.seed(1)
    n_member <- 2
    n_time <- 4
    n_group <- 30
    dd <- make_sep_dat(n_member = n_member, n_time = n_time, n_group = n_group)
    sd <- c(0.8, 1.2)
    rho <- 0.25
    phi <- 0.4
    sigma <- 0.5
    R_member <- cs_corr(n_member, rho)
    R_time <- ar1_corr(n_time, phi)
    sd_full <- rep(sd, n_time)
    Sigma <- diag(sd_full) %*% kronecker(R_time, R_member) %*% diag(sd_full)
    B <- t(matrix(rnorm(n_group * n_member * n_time), nrow = n_group) %*%
               chol(Sigma))
    dd$y <- as.vector(B) + rnorm(nrow(dd), sd = sigma)

    fit <- glmmTMB(y ~ 1 +
                       separable(us(0 + member) %x% ar1(0 + time) | group),
                   data = dd)

    expect_equal(fit$fit$convergence, 0)
})

test_that("separable prediction with newdata reports current limitation", {
    dd <- make_sep_dat()
    theta <- c(log(1), qlogis((0.2 + 1) / 2), ar1_to_theta(0.3))
    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    expect_error(predict(fit, newdata = dd[1, ]),
                 "newdata is not yet implemented")
})

test_that("separable simulation works for product covariance structures", {
    dd <- make_sep_dat(n_time = 2, reps = TRUE)
    theta <- c(log(1), qlogis((0.2 + 1) / 2), ar1_to_theta(0.3))
    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    sims <- simulate(fit, nsim = 2)
    expect_s3_class(sims, "data.frame")
    expect_equal(dim(sims), c(nrow(dd), 2L))
    expect_true(all(vapply(sims, is.numeric, logical(1))))
})
