stopifnot(require("testthat"), require("glmmTMB"))

## Contract tests for the minimal Kronecker API.  The intended parameter order is
## one global log-SD followed by the normalized margin parameters, left to right.

kron_test_data <- function(n_group = 2L) {
    dd <- expand.grid(
        member = factor(paste0("m", 1:2)),
        time = factor(paste0("t", 1:3)),
        group = factor(seq_len(n_group))
    )
    dd$y <- 0
    dd
}

kron_row_product <- function(...) {
    X <- list(...)
    ans <- X[[1L]]
    if (length(X) > 1L) {
        for (i in seq.int(2L, length(X))) {
            ans <- do.call(cbind, lapply(seq_len(ncol(X[[i]])), function(j) {
                ans * X[[i]][, j]
            }))
        }
    }
    colnames(ans) <- Reduce(
        function(a, b) as.vector(outer(a, b, paste, sep = ":")),
        lapply(X, colnames)
    )
    ans
}

kron_slope_design <- function(data) {
    kron_row_product(
        model.matrix(~ 0 + member + member:x, data),
        model.matrix(~ 0 + time, data)
    )
}

kron_ar1_theta <- function(phi) phi / sqrt(1 - phi^2)

kron_homcs_theta <- function(rho, n) {
    lower <- -1 / (n - 1)
    qlogis((rho - lower) / (1 - lower))
}

kron_ar1_phi <- function(theta) theta / sqrt(1 + theta^2)

kron_ar1_cor <- function(n, phi) {
    outer(seq_len(n), seq_len(n), function(i, j) phi^abs(i - j))
}

kron_cs_cor <- function(n, rho) {
    ans <- matrix(rho, n, n)
    diag(ans) <- 1
    ans
}

kron_fix_theta <- function(formula, data, theta, dispformula = ~1) {
    glmmTMB(
        formula, data = data, dispformula = dispformula,
        start = list(theta = theta),
        map = list(theta = factor(rep(NA, length(theta)))),
        control = glmmTMBControl(conv_check = "skip",
                                 eigval_check = FALSE)
    )
}

kron_joint_state <- function(formula, data, theta, b, residual_sd = 2) {
    fit <- kron_fix_theta(formula, data, theta)
    parameters <- fit$obj$env$parList(fit$fit$par, fit$fit$parfull)
    joint <- TMB::MakeADFun(fit$obj$env$data, parameters,
                            random = NULL, DLL = "glmmTMB", silent = TRUE)
    par <- joint$par
    par[names(par) == "beta"] <- 0
    par[names(par) == "betadisp"] <- log(residual_sd)
    par[names(par) == "b"] <- b
    par[names(par) == "theta"] <- theta
    list(fit = fit, obj = joint, par = par)
}

kron_joint_eval <- function(state, theta = NULL, b = NULL) {
    par <- state$par
    if (!is.null(theta)) par[names(par) == "theta"] <- theta
    if (!is.null(b)) par[names(par) == "b"] <- b
    state$obj$fn(par)
}

test_that("kron builds the exact random-slope product design", {
    dd <- kron_test_data()
    dd$x <- seq(-1, 1, length.out = nrow(dd))
    dd <- dd[!(dd$member == "m2" & dd$time == "t2"), ]
    prep <- glmmTMB(
        y ~ 1 + kron(us(0 + member + member:x) %x%
                         ar1(0 + time) | group),
        data = dd, doFit = FALSE
    )

    product <- kron_slope_design(dd)
    expected <- do.call(rbind, lapply(levels(dd$group), function(g) {
        t(product * (dd$group == g))
    }))
    expect_equal(
        as.matrix(prep$condList$reTrms$Zt), as.matrix(expected),
        check.attributes = FALSE
    )
    expect_equal(prep$condReStruc[[1L]]$kronDims, c(4L, 3L))
    expect_equal(prep$condReStruc[[1L]]$kronMarginColumns,
                 list(colnames(model.matrix(~ 0 + member + member:x, dd)),
                      colnames(model.matrix(~ 0 + time, dd))))
    expect_equal(prep$condList$reTrms$cnms[[1L]], colnames(product))
    expect_length(prep$parameters$theta, 11L)

    conventional <- glmmTMB(
        y ~ 1 + kron(us(1 + x) %x% ar1(0 + time) | group),
        data = dd, doFit = FALSE
    )
    expect_equal(conventional$condReStruc[[1L]]$kronDims, c(2L, 3L))
    expect_equal(conventional$condReStruc[[1L]]$kronMarginColumns[[1L]],
                 c("(Intercept)", "x"))

    dd$period <- rep(1:3, length.out = nrow(dd))
    period_levels <- 1:5
    inline_factor <- glmmTMB(
        y ~ kron(homdiag(0 + factor(period, levels = period_levels)) %x%
                     ar1(0 + time) | group),
        data = dd, doFit = FALSE
    )
    expect_equal(inline_factor$condReStruc[[1L]]$kronDims, c(5L, 3L))
})

test_that("kron random slopes match dense us likelihood and gradients", {
    skip_if_not_installed("numDeriv")
    dd <- kron_test_data()
    dd$x <- seq(-0.8, 0.9, length.out = nrow(dd))
    dd$P <- I(kron_slope_design(dd))

    phi <- 0.45
    global_sd <- 1.3
    logsd <- c(0.2, -0.1, 0.05)
    relative_sd <- exp(c(logsd, -sum(logsd)))
    us_raw <- c(0.15, -0.1, 0.05, 0.12, -0.08, 0.1)
    R_slope <- get_cor(us_raw, "mat")
    Sigma_slope <- outer(relative_sd, relative_sd) * R_slope
    R_time <- kron_ar1_cor(3, phi)
    Sigma_full <- global_sd^2 * kronecker(R_time, Sigma_slope)
    theta <- c(log(global_sd), logsd, us_raw, kron_ar1_theta(phi))
    dense_theta <- c(log(sqrt(diag(Sigma_full))),
                     put_cor(cov2cor(Sigma_full)))
    b <- seq(-0.4, 0.5, length.out = 24)

    kron_form <- y ~ 1 + kron(us(0 + member + member:x) %x%
                                  ar1(0 + time) | group)
    dense_form <- y ~ 1 + us(0 + P | group)

    ks <- kron_joint_state(kron_form, dd, theta, b)
    ds <- kron_joint_state(dense_form, dd, dense_theta, b)

    dense_from_kron <- function(x) {
        sd1 <- exp(c(x[2:4], -sum(x[2:4])))
        Sigma1 <- outer(sd1, sd1) * get_cor(x[5:10], "mat")
        Sigma2 <- kron_ar1_cor(3, kron_ar1_phi(x[11]))
        Sigma <- exp(2 * x[1]) * kronecker(Sigma2, Sigma1)
        c(log(sqrt(diag(Sigma))), put_cor(cov2cor(Sigma)))
    }
    fk <- function(x) kron_joint_eval(ks, theta = x)
    fd <- function(x) kron_joint_eval(ds, theta = dense_from_kron(x))

    expect_equal(unname(fk(theta)), unname(fd(theta)), tolerance = 1e-6)
    expect_equal(numDeriv::grad(fk, theta), numDeriv::grad(fd, theta),
                 tolerance = 2e-5)

    fbk <- function(x) kron_joint_eval(ks, b = x)
    fbd <- function(x) kron_joint_eval(ds, b = x)
    expect_equal(numDeriv::grad(fbk, b), numDeriv::grad(fbd, b),
                 tolerance = 2e-5)

    vc <- VarCorr(ks$fit, full_cor = TRUE)$cond[[1]]
    expect_equal(unname(attr(vc, "stddev")), sqrt(diag(Sigma_full)),
                 tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")), cov2cor(Sigma_full),
                 tolerance = 1e-6)

    prep <- glmmTMB(kron_form, data = dd, doFit = FALSE)
    expect_length(prep$parameters$theta, 11L)
})

test_that("kron accepts more than five margins", {
    dd <- expand.grid(
        u1 = factor(1:2), u2 = factor(1:2), u3 = factor(1:2),
        u4 = factor(1:2), u5 = factor(1:2), u6 = factor(1:2),
        group = factor(1:2)
    )
    set.seed(101)
    dd$y <- rnorm(nrow(dd))
    theta <- c(log(0.8), rep(kron_ar1_theta(0.15), 6))

    fit <- kron_fix_theta(
        y ~ 1 + kron(ar1(0 + u1) %x% ar1(0 + u2) %x%
                         ar1(0 + u3) %x% ar1(0 + u4) %x%
                         ar1(0 + u5) %x% ar1(0 + u6) | group),
        dd, theta
    )

    expect_true(is.finite(as.numeric(logLik(fit))))
    expect_equal(length(fit$obj$env$parList(
        fit$fit$par, fit$obj$env$last.par.best)$theta), 7L)
})

test_that("heterogeneous margins use identifiable geometric-mean normalization", {
    skip_if_not_installed("numDeriv")
    dd <- expand.grid(
        axis1 = factor(1:3), axis2 = factor(1:2), group = factor(1:4)
    )
    dd$y <- 0
    theta <- c(log(1.2), 0.25, -0.15, 0.2)
    b <- seq(-0.7, 0.8, length.out = 24)
    form <- y ~ 1 + kron(diag(0 + axis1) %x%
                             diag(0 + axis2) | group)
    state <- kron_joint_state(form, dd, theta, b)

    ## One global log-SD plus d-1 contrasts for each diagonal margin.
    prep <- glmmTMB(form, data = dd, doFit = FALSE)
    expect_length(prep$parameters$theta, 4L)

    cell_sd <- unname(attr(VarCorr(state$fit)$cond[[1]], "stddev"))
    relative_sd <- cell_sd / exp(theta[1])
    expect_equal(exp(mean(log(relative_sd))), 1, tolerance = 1e-7)

    ## WHITE-BOX/NUMERICAL: the old unconstrained product parameterization had
    ## an exact scale ridge.  The four-parameter Hessian must now have full rank.
    H <- numDeriv::hessian(
        function(x) kron_joint_eval(state, theta = x), theta
    )
    expect_equal(qr(H, tol = 1e-5)$rank, 4L)
})

test_that("every supported margin matches a dense covariance", {
    dd <- expand.grid(
        axis1 = factor(1:3), axis2 = factor(1:2), group = factor(1:2)
    )
    dd$y <- 0
    dd$P <- I(kron_row_product(
        model.matrix(~ 0 + axis1, dd),
        model.matrix(~ 0 + axis2, dd)
    ))
    b <- seq(-0.5, 0.6, length.out = 12)
    global_sd <- 1.15
    logsd <- c(0.2, -0.1)
    relative_sd <- exp(c(logsd, -sum(logsd)))
    rho <- 0.15
    phi <- 0.3
    us_raw <- c(0.2, -0.15, 0.1)

    with_sd <- function(corr) outer(relative_sd, relative_sd) * corr
    cases <- list(
        homdiag = list(theta = numeric(), covariance = diag(3)),
        diag = list(theta = logsd, covariance = diag(relative_sd^2)),
        homcs = list(theta = kron_homcs_theta(rho, 3),
                     covariance = kron_cs_cor(3, rho)),
        cs = list(theta = c(logsd, kron_homcs_theta(rho, 3)),
                  covariance = with_sd(kron_cs_cor(3, rho))),
        us = list(theta = c(logsd, us_raw),
                  covariance = with_sd(get_cor(us_raw, "mat"))),
        ar1 = list(theta = kron_ar1_theta(phi),
                   covariance = kron_ar1_cor(3, phi)),
        hetar1 = list(theta = c(logsd, kron_ar1_theta(phi)),
                      covariance = with_sd(kron_ar1_cor(3, phi)))
    )

    for (margin in names(cases)) {
        form <- as.formula(sprintf(
            "y ~ 1 + kron(%s(0 + axis1) %%x%% homdiag(0 + axis2) | group)",
            margin
        ))
        theta <- c(log(global_sd), cases[[margin]]$theta)
        full_cov <- global_sd^2 * kronecker(diag(2),
                                             cases[[margin]]$covariance)
        dense_theta <- c(log(sqrt(diag(full_cov))),
                         put_cor(cov2cor(full_cov)))
        ks <- kron_joint_state(form, dd, theta, b)
        ds <- kron_joint_state(
            y ~ 1 + us(0 + P | group),
            dd, dense_theta, b
        )

        expect_equal(kron_joint_eval(ks), kron_joint_eval(ds),
                     tolerance = 1e-6, info = margin)
        vc <- VarCorr(ks$fit, full_cor = TRUE)$cond[[1]]
        vc_values <- matrix(as.numeric(vc), nrow = nrow(full_cov))
        expect_equal(vc_values, unname(full_cov), tolerance = 1e-6,
                     info = margin)
    }

    expect_error(
        glmmTMB(
            y ~ kron(toep(0 + axis1) %x%
                         homdiag(0 + axis2) | group),
            dd, doFit = FALSE
        ),
        "unsupported kron.*toep"
    )
})

test_that("AR1 product precision remains sparse", {
    dd <- expand.grid(a = factor(1:2), time = factor(1:8),
                      group = factor(1))
    dd$y <- 0
    prep <- glmmTMB(
        y ~ 1 + kron(homdiag(0 + a) %x%
                         ar1(0 + time) | group),
        dd, start = list(theta = c(0, 0.5)), doFit = FALSE
    )
    obj <- fitTMB(prep, doOptim = FALSE)
    random_hessian <- obj$env$spHess(obj$env$par, random = TRUE)

    ## Two independent AR1(8) precision blocks: 8 diagonal + 7 off-diagonal
    ## entries per upper triangle. A dense marginal would store 72 entries.
    expect_equal(dim(random_hessian), c(16L, 16L))
    expect_equal(length(random_hessian@x), 2L * (8L + 7L))
})

test_that("kron simulation has the requested product covariance", {
    dd <- kron_test_data(n_group = 1L)
    rho <- 0.25
    phi <- 0.4
    global_sd <- 1.1
    theta <- c(log(global_sd),
               kron_homcs_theta(rho, 2), kron_ar1_theta(phi))
    fit <- kron_fix_theta(
        y ~ 1 + kron(homcs(0 + member) %x%
                         ar1(0 + time) | group),
        dd, theta, dispformula = ~0
    )

    sim1 <- simulate(fit, nsim = 1000, seed = 2026)
    sim2 <- simulate(fit, nsim = 1000, seed = 2026)
    expect_equal(sim1, sim2)

    expected <- global_sd^2 * kronecker(
        kron_ar1_cor(3, phi), kron_cs_cor(2, rho)
    )
    empirical <- cov(t(as.matrix(sim1)))
    expect_lt(max(abs(empirical - expected)), 0.15)
})

test_that("kron covariance parameters can be freely optimized", {
    set.seed(901)
    dd <- kron_test_data(n_group = 20L)
    covariance <- 1.2^2 * kronecker(
        kron_ar1_cor(3, 0.45), kron_cs_cor(2, 0.25)
    )
    random_effect <- as.vector(
        t(chol(covariance)) %*% matrix(rnorm(120), 6L)
    )
    dd$y <- 1 + random_effect + rnorm(nrow(dd), sd = 0.5)

    fit <- glmmTMB(
        y ~ 1 + kron(homcs(0 + member) %x%
                         ar1(0 + time) | group),
        dd, start = list(theta = c(0, 0, 0))
    )

    theta <- fit$obj$env$parList(fit$fit$par, fit$fit$parfull)$theta
    expect_equal(fit$fit$convergence, 0)
    expect_true(fit$sdr$pdHess)
    expect_true(all(is.finite(theta)))
    expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("kron random slopes remain usable by update and predict", {
    dd <- kron_test_data(n_group = 6L)
    dd$x <- seq(-1, 1, length.out = nrow(dd))
    dd$y <- 1 + 0.5 * dd$x
    theta <- c(log(0.7), rep(0, 9), kron_ar1_theta(0.3))
    user_formula <- y ~ 1 +
        kron(us(0 + member + member:x) %x% ar1(0 + time) | group)
    ## Call glmmTMB directly here: update() must exercise the stored user call,
    ## not a helper call whose `data` promise has gone out of scope.
    fit <- glmmTMB(
        user_formula, data = dd,
        start = list(theta = theta),
        map = list(theta = factor(rep(NA, length(theta)))),
        control = glmmTMBControl(conv_check = "skip")
    )

    expect_equal(formula(fit), user_formula)
    expect_equal(formula(fit, component = "cond"), user_formula)

    updated <- update(fit, . ~ . + x)
    form_text <- paste(deparse(formula(updated)), collapse = " ")
    expect_match(form_text, "kron\\(")
    expect_false(grepl("\\.\\.(kron|sep)|data.frame", form_text))

    train_pop <- predict(fit, re.form = NA)
    expect_equal(predict(fit, newdata = dd, re.form = NA), train_pop)
    pop_without_re <- predict(
        fit, newdata = data.frame(dummy = 1:3), re.form = NA
    )
    expect_equal(pop_without_re, rep(train_pop[[1L]], 3L))
    train_cond <- predict(fit)
    expect_gt(diff(range(train_cond)), 0.1)

    take <- c(1L, 8L, 15L)
    newdata <- droplevels(dd[take, ])
    expect_equal(predict(fit, newdata = newdata), train_cond[take])

    changed_x <- newdata
    changed_x$x <- changed_x$x + 0.75
    expect_gt(max(abs(predict(fit, newdata = changed_x) -
                      train_cond[take])), 1e-4)

    ## Margin labels, not the integer codes of newdata factors, identify columns.
    reversed <- dd
    reversed$member <- factor(as.character(reversed$member),
                              levels = rev(levels(dd$member)))
    reversed$time <- factor(as.character(reversed$time),
                            levels = rev(levels(dd$time)))
    expect_equal(predict(fit, newdata = reversed), train_cond)

    unseen <- dd[1, ]
    unseen$member <- factor("new-member")
    expect_error(predict(fit, newdata = unseen, allow.new.levels = TRUE),
                 "fitted|new level|unknown")

    new_group <- dd[1:2, ]
    new_group$group <- factor("new-group")
    expect_true(all(is.finite(predict(
        fit, newdata = new_group, allow.new.levels = TRUE
    ))))

    isolated_formula <- as.formula(
        paste0("y ~ kron(us(0 + member + member:x) %x% ",
               "ar1(0 + time) | group)"),
        env = new.env(parent = baseenv())
    )
    expect_no_error(glmmTMB(isolated_formula, data = dd, doFit = FALSE))
})

test_that("kron prediction reuses fitted marginal bases", {
    dd <- kron_test_data(n_group = 4L)
    dd$x <- seq(-1, 1, length.out = nrow(dd))
    dd$y <- 1 + 0.5 * dd$x + 0.8 * dd$x^2
    degree <- 2
    theta <- c(log(0.7), 0, 0, kron_ar1_theta(0.2))
    fit <- glmmTMB(
        y ~ 1 + kron(us(0 + poly(x, degree)) %x%
                         ar1(0 + time) | group),
        data = dd, start = list(theta = theta),
        map = list(theta = factor(rep(NA, length(theta)))),
        control = glmmTMBControl(conv_check = "skip")
    )

    train <- predict(fit)
    augmented <- rbind(dd, transform(dd[1L, ], x = 50))
    pred <- predict(fit, newdata = augmented)
    expect_equal(pred[seq_len(nrow(dd))], train, tolerance = 1e-10)
    expect_true(is.finite(tail(pred, 1L)))
})

test_that("population prediction works after fitting without data", {
    dd <- kron_test_data(n_group = 3L)
    x <- seq(-1, 1, length.out = nrow(dd))
    time <- dd$time
    group <- dd$group
    y <- 1 + x
    degree <- 2
    theta <- c(log(0.7), 0, 0, kron_ar1_theta(0.2))
    fit <- suppressWarnings(glmmTMB(
        y ~ 1 + kron(us(0 + poly(x, degree)) %x%
                         ar1(0 + time) | group),
        start = list(theta = theta),
        map = list(theta = factor(rep(NA, length(theta)))),
        control = glmmTMBControl(conv_check = "skip")
    ))

    expect_equal(fit$modelInfo$reStruc$condReStruc[[1L]]$kronSourceVars,
                 c("x", "time"))
    expect_true(all(is.finite(predict(
        fit, newdata = data.frame(row = 1:2), re.form = NA
    ))))
})

test_that("kron terms compose across model terms and components", {
    dd <- kron_test_data()
    dd$x <- seq(-1, 1, length.out = nrow(dd))
    dd$batch <- factor(rep(1:2, length.out = nrow(dd)))

    prep <- glmmTMB(
        y ~ 1 + (1 | batch) +
            kron(homcs(0 + member) %x% ar1(0 + time) | group) +
            kron(us(0 + member + member:x) %x%
                     homdiag(0 + time) | batch),
        ziformula = ~kron(us(0 + member + member:x) %x%
                              ar1(0 + time) | group),
        dispformula = ~kron(homdiag(0 + member) %x%
                                diag(0 + time) | batch),
        data = dd, doFit = FALSE
    )

    expect_equal(prep$condList$ss, c("us", "kron", "kron"))
    expect_equal(unname(lapply(prep$condReStruc[-1L], `[[`, "kronDims")),
                 list(c(2L, 3L), c(4L, 3L)))
    expect_equal(prep$ziReStruc[[1L]]$kronDims, c(4L, 3L))
    expect_equal(prep$dispReStruc[[1L]]$kronDims, c(2L, 3L))
    expect_equal(unname(prep$condList$reTrms$Gp),
                 c(0L, 2L, 14L, 38L))

    obj <- fitTMB(prep, doOptim = FALSE)
    expect_true(is.finite(obj$fn(obj$par)))
    expect_true(all(is.finite(obj$gr(obj$par))))
})

test_that("fixed-effect Kronecker products are left alone", {
    dd <- kron_test_data()
    A <- matrix(seq_len(nrow(dd)), ncol = 1L)
    prep <- glmmTMB(
        y ~ I(A %x% 1) +
            kron(homdiag(0 + member) %x% ar1(0 + time) | group),
        data = dd, doFit = FALSE
    )

    expect_equal(ncol(prep$condList$X), 2L)
    expect_equal(unname(prep$condList$X[, 2L]), as.vector(A))
})

test_that("Kronecker reporting is factorized by default", {
    dd <- kron_test_data()
    theta <- c(log(0.8), kron_homcs_theta(0.1, 2),
               kron_ar1_theta(0.2))
    fit <- glmmTMB(
        y ~ 1 + kron(homcs(0 + member) %x%
                         ar1(0 + time) | group),
        data = dd, start = list(theta = theta),
        map = list(theta = factor(rep(NA, length(theta))))
    )
    vc <- VarCorr(fit)$cond[[1]]
    factors <- attr(vc, "kron")

    expect_true(is.na(vc))
    expect_equal(factors$dims, c(2L, 3L))
    expect_equal(dim(factors$margins[[1L]]), c(2L, 2L))
    expect_equal(dim(factors$margins[[2L]]), c(3L, 3L))
    expect_false(any(vapply(factors$margins, nrow, integer(1)) == 6L))

    full <- VarCorr(fit, full_cor = TRUE)$cond[[1]]
    expect_equal(dim(full), c(6L, 6L))
})
