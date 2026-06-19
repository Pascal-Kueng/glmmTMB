stopifnot(require("testthat"),
          require("glmmTMB"))

make_sep_dat <- function(n_member = 2, n_time = 3, n_group = 2, reps = FALSE) {
    args <- list(member = factor(paste0("m", seq_len(n_member))),
                 time = factor(seq_len(n_time)),
                 group = factor(seq_len(n_group)))
    if (reps) args$rep <- 1:2
    dd <- do.call(expand.grid, args)
    dd$y <- seq_len(nrow(dd)) / nrow(dd)
    dd
}

ar1_to_theta <- function(phi) phi / sqrt(1 - phi^2)

homcs_to_theta <- function(rho, n_member) {
    lower <- -1 / (n_member - 1)
    qlogis((rho - lower) / (1 - lower))
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

fit_fixed_theta <- function(form, dd, theta) {
    glmmTMB(form, data = dd,
            start = list(theta = theta),
            map = list(theta = factor(rep(NA, length(theta)))))
}

make_sep_case <- function(struc = c("cs", "homcs", "us"), reversed = FALSE,
                          n_member = 2, n_time = 3,
                          scale_mode = c("margin", "global", "product",
                                         "selected_product")) {
    struc <- match.arg(struc)
    scale_mode <- match.arg(scale_mode)
    rho <- switch(struc, cs = 0.2, homcs = 0.3, us = -0.25)
    phi <- switch(struc, cs = 0.45, homcs = 0.4, us = 0.5)
    sd <- switch(struc,
        cs = seq(0.8, 1.2, length.out = n_member),
        homcs = 2,
        us = seq(0.8, 1.2, length.out = n_member)
    )
    global_sd <- 1.4

    R_member <- if (struc %in% c("cs", "homcs")) {
        M <- matrix(rho, n_member, n_member)
        diag(M) <- 1
        M
    } else {
        outer(seq_len(n_member), seq_len(n_member),
              function(i, j) ifelse(i == j, 1, 0.35^abs(i - j)))
    }
    R_time <- outer(seq_len(n_time), seq_len(n_time),
                    function(i, j) phi^abs(i - j))
    member_scale_theta <- switch(struc,
        cs = log(sd),
        homcs = log(sd),
        us = log(sd)
    )
    member_corr_theta <- switch(struc,
        cs = homcs_to_theta(rho, n_member),
        homcs = homcs_to_theta(rho, n_member),
        us = put_cor(R_member)
    )
    member_theta <- if (scale_mode == "global") {
        member_corr_theta
    } else {
        c(member_scale_theta, member_corr_theta)
    }
    scale_call <- switch(scale_mode,
        margin = NULL,
        global = quote(global()),
        product = quote(product()),
        selected_product = as.call(list(as.name("product"),
                                     as.call(list(as.name(struc),
                                                  quote(0 + member)))))
    )
    sep_call <- function(lhs, group = quote(group)) {
        args <- list(as.name("separable"),
                     as.call(list(as.name("|"), lhs, group)))
        if (!is.null(scale_call)) args$scale <- scale_call
        as.call(args)
    }
    dense_call <- as.call(list(as.name(struc), quote(0 + member)))
    ar1_call <- quote(ar1(0 + time))

    if (reversed) {
        form <- as.formula(as.call(list(quote(`~`), quote(y),
            as.call(list(quote(`+`), 1,
                         sep_call(as.call(list(as.name("%x%"),
                                               ar1_call, dense_call))))))))
        dense_form <- y ~ 1 + us(sepgrid(time, member) + 0 | group)
        theta <- if (scale_mode == "global") {
            c(log(global_sd), ar1_to_theta(phi), member_theta)
        } else {
            c(ar1_to_theta(phi), member_theta)
        }
        R_full <- kronecker(R_member, R_time)
        sd_full <- if (scale_mode == "global") {
            rep(global_sd, n_member * n_time)
        } else if (length(sd) == 1) {
            rep(sd, n_member * n_time)
        } else {
            rep(sd, each = n_time)
        }
        codes <- unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[[struc]]))
        kinds <- c(2L, 1L)
        scale_spec <- if (scale_mode == "global") integer() else 1L
    } else {
        form <- as.formula(as.call(list(quote(`~`), quote(y),
            as.call(list(quote(`+`), 1,
                         sep_call(as.call(list(as.name("%x%"),
                                               dense_call, ar1_call))))))))
        dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
        theta <- if (scale_mode == "global") {
            c(log(global_sd), member_theta, ar1_to_theta(phi))
        } else {
            c(member_theta, ar1_to_theta(phi))
        }
        R_full <- kronecker(R_time, R_member)
        sd_full <- if (scale_mode == "global") {
            rep(global_sd, n_member * n_time)
        } else if (length(sd) == 1) {
            rep(sd, n_member * n_time)
        } else {
            rep(sd, n_time)
        }
        codes <- unname(c(.valid_covstruct[[struc]], .valid_covstruct[["ar1"]]))
        kinds <- c(1L, 2L)
        scale_spec <- if (scale_mode == "global") integer() else 0L
    }

    scale_mode_code <- switch(scale_mode,
        margin = 1L, global = 2L, product = 3L, selected_product = 4L)
    list(form = form, dense_form = dense_form, theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full, sd_full = sd_full,
         codes = codes, kinds = kinds, dispatch = 1L,
         scale_mode = scale_mode_code,
         scale_spec = scale_spec,
         n_member = n_member, n_time = n_time)
}

case_data <- function(case, reps = FALSE, n_group = 2) {
    if (!is.null(case$dd)) return(case$dd)
    make_sep_dat(n_member = case$n_member, n_time = case$n_time,
                 n_group = n_group, reps = reps)
}

expect_separable_case_vc <- function(case) {
    dd <- case_data(case, reps = TRUE)
    fit <- fit_fixed_theta(case$form, dd, case$theta)
    restruc <- fit$modelInfo$reStruc$condReStruc[[1]]
    vc <- VarCorr(fit)$cond[[1]]

    expect_equal(restruc$sepCodes, case$codes)
    expect_equal(restruc$sepDensityKinds, case$kinds)
    expect_equal(restruc$sepDispatch, case$dispatch)
    expect_equal(restruc$sepScaleMode, case$scale_mode)
    expect_equal(restruc$sepScaleSpec, case$scale_spec)
    expect_equal(unname(attr(vc, "stddev")), case$sd_full, tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")), case$R_full, tolerance = 1e-6)
}

expect_separable_case_nll <- function(case) {
    dd <- case_data(case, n_group = 2)
    dd$y <- 0
    b <- seq(-0.4, 0.5, length.out = length(case$sd_full) * nlevels(dd$group))

    sep_nll <- joint_nll_at(case$form, dd, case$theta, b)
    dense_nll <- joint_nll_at(case$dense_form, dd, case$theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
}

make_dense_margin <- function(struc = c("cs", "homcs", "us"), n = 2,
                              sd = seq(0.8, 1.2, length.out = n),
                              rho = 0.25) {
    struc <- match.arg(struc)
    R <- if (struc %in% c("cs", "homcs")) {
        M <- matrix(rho, n, n)
        diag(M) <- 1
        M
    } else {
        outer(seq_len(n), seq_len(n),
              function(i, j) ifelse(i == j, 1, rho^abs(i - j)))
    }
    list(
        struc = struc,
        sd = if (struc == "homcs") sd[[1]] else sd,
        R = R,
        scale_theta = log(if (struc == "homcs") sd[[1]] else sd),
        corr_theta = if (struc %in% c("cs", "homcs")) homcs_to_theta(rho, n)
                     else put_cor(R)
    )
}

make_sep_dense_dense_case <- function(struc0 = c("cs", "homcs", "us"),
                                      struc1 = c("cs", "homcs", "us"),
                                      scale_mode = c("global", "product",
                                                     "selected_first",
                                                     "selected_second"),
                                      n0 = 2, n1 = 3) {
    struc0 <- match.arg(struc0)
    struc1 <- match.arg(struc1)
    scale_mode <- match.arg(scale_mode)
    m0 <- make_dense_margin(struc0, n0, rho = 0.2)
    m1 <- make_dense_margin(struc1, n1, sd = seq(1.1, 1.4, length.out = n1),
                            rho = 0.35)
    global_sd <- 1.3

    dd <- expand.grid(member = factor(paste0("m", seq_len(n0))),
                      item = factor(paste0("i", seq_len(n1))),
                      group = factor(seq_len(2)))
    dd$y <- 0

    margin_call <- function(struc, var) {
        as.call(list(as.name(struc),
                     as.call(list(as.name("+"), 0, as.name(var)))))
    }
    call0 <- margin_call(struc0, "member")
    call1 <- margin_call(struc1, "item")
    scale_call <- switch(scale_mode,
        global = quote(global()),
        product = quote(product()),
        selected_first = as.call(list(as.name("product"), call0)),
        selected_second = as.call(list(as.name("product"), call1))
    )
    sep_call <- as.call(list(
        as.name("separable"),
        as.call(list(as.name("|"),
                     as.call(list(as.name("%x%"), call0, call1)),
                     quote(group))),
        scale = scale_call
    ))
    form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1, sep_call)))))

    theta0 <- switch(scale_mode,
        global = m0$corr_theta,
        product = c(m0$scale_theta, m0$corr_theta),
        selected_first = c(m0$scale_theta, m0$corr_theta),
        selected_second = m0$corr_theta
    )
    theta1 <- switch(scale_mode,
        global = m1$corr_theta,
        product = c(m1$scale_theta, m1$corr_theta),
        selected_first = m1$corr_theta,
        selected_second = c(m1$scale_theta, m1$corr_theta)
    )
    theta <- if (scale_mode == "global") {
        c(log(global_sd), theta0, theta1)
    } else {
        c(theta0, theta1)
    }

    sd0 <- if (length(m0$sd) == 1L) rep(m0$sd, n0) else m0$sd
    sd1 <- if (length(m1$sd) == 1L) rep(m1$sd, n1) else m1$sd
    sd_full <- switch(scale_mode,
        global = rep(global_sd, n0 * n1),
        product = as.vector(outer(sd0, sd1)),
        selected_first = rep(sd0, n1),
        selected_second = rep(sd1, each = n0)
    )
    R_full <- kronecker(m1$R, m0$R)

    list(form = form,
         dense_form = y ~ 1 + us(sepgrid(member, item) + 0 | group),
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = unname(c(.valid_covstruct[[struc0]],
                          .valid_covstruct[[struc1]])),
         kinds = c(1L, 1L),
         dispatch = 2L,
         scale_mode = switch(scale_mode,
             global = 2L, product = 3L,
             selected_first = 4L, selected_second = 4L),
         scale_spec = switch(scale_mode,
             global = integer(), product = c(0L, 1L),
             selected_first = 0L, selected_second = 1L))
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

test_that("glmmTMB preserves unused separable margin levels by default", {
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

test_that("separable specs handle product order", {
    dd <- make_sep_dat()

    h <- glmmTMB(y ~ 1 +
                     separable(ar1(0 + time) %x% homcs(0 + member) | group),
                 data = dd, doFit = FALSE)
    u <- glmmTMB(y ~ 1 +
                     separable(us(0 + member) %x% ar1(0 + time) | group),
                 data = dd, doFit = FALSE)

    expect_equal(unname(h$condReStruc[[1]]$blockCode),
                 unname(.valid_covstruct[["separable"]]))
    expect_equal(h$condReStruc[[1]]$blockNumTheta, 3)
    expect_equal(h$condReStruc[[1]]$sepDims, c(3L, 2L))
    expect_equal(h$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[["homcs"]])))
    expect_equal(h$condReStruc[[1]]$sepDensityKinds, c(2L, 1L))
    expect_equal(h$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleSpec, 1L)

    expect_equal(u$condReStruc[[1]]$blockNumTheta, 4)
    expect_equal(u$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["us"]], .valid_covstruct[["ar1"]])))
    expect_equal(u$condReStruc[[1]]$sepDensityKinds, c(1L, 2L))
    expect_equal(u$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(u$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(u$condReStruc[[1]]$sepScaleSpec, 0L)
})

test_that("separable parser stores structured specs outside splitForm payload", {
    f <- y ~ 1 +
        separable(homcs(0 + member) %x% ar1(0 + time) | group,
                  scale = homcs(0 + member))
    g <- glmmTMB:::rewrite_separable_formula(f)
    ss <- reformulas::splitForm(g, specials = c(names(.valid_covstruct), "s"))
    specs <- attr(g, "separable_specs")
    spec <- specs[[1]]

    expect_length(specs, 1)
    expect_equal(spec$grid, c("member", "time"))
    expect_equal(unname(spec$margins$struc), c("homcs", "ar1"))
    expect_equal(spec$scale$mode, "margin")
    expect_equal(unname(spec$scale$margins$struc), "homcs")
    expect_equal(unname(spec$scale$margins$var), "member")
    expect_equal(ss$reTrmClasses, "separable")
    expect_equal(deparse(ss$reTrmFormulas[[1]]),
                 "0 + (0 + member + (0 + time)) | group")
    expect_equal(length(ss$reTrmAddArgs[[1]]), 2)
    expect_equal(eval(ss$reTrmAddArgs[[1]][[2]]), 1L)
})

test_that("separable parser flattens product chains and records scale syntax", {
    f <- y ~ 1 +
        separable(us(0 + member) %x% ar1(0 + time) %x% cs(0 + item) | group,
                  scale = product(us(0 + member), cs(0 + item)))
    g <- glmmTMB:::rewrite_separable_formula(f)
    spec <- attr(g, "separable_specs")[[1]]

    expect_equal(spec$grid, c("member", "time", "item"))
    expect_equal(spec$margins$struc, c("us", "ar1", "cs"))
    expect_equal(spec$scale$mode, "selected_product")
    expect_equal(spec$scale$margins$struc, c("us", "cs"))
    expect_equal(deparse(g[[3]][[3]][[2]]),
                 "0 + (0 + member + (0 + time) + (0 + item)) | group")
})

test_that("separable frontend parses simple existing covariance margins", {
    f <- y ~ 1 +
        separable(diag(0 + member) %x% ar1(0 + time) %x%
                      homtoep(0 + item) | group,
                  scale = product(diag(0 + member), homtoep(0 + item)))
    g <- glmmTMB:::rewrite_separable_formula(f)
    spec <- attr(g, "separable_specs")[[1]]

    expect_equal(spec$grid, c("member", "time", "item"))
    expect_equal(spec$margins$struc, c("diag", "ar1", "homtoep"))
    expect_equal(spec$scale$mode, "selected_product")
    expect_equal(spec$scale$margins$struc, c("diag", "homtoep"))
})

test_that("separable parser records global and product scale modes", {
    f_global <- y ~ 1 +
        separable(ar1(0 + member) %x% ar1(0 + time) | group,
                  scale = global())
    g_global <- glmmTMB:::rewrite_separable_formula(f_global)
    expect_equal(attr(g_global, "separable_specs")[[1]]$scale$mode, "global")

    f_product <- y ~ 1 +
        separable(us(0 + member) %x% cs(0 + item) | group,
                  scale = product())
    g_product <- glmmTMB:::rewrite_separable_formula(f_product)
    expect_equal(attr(g_product, "separable_specs")[[1]]$scale$mode, "product")
})

test_that("separable selected product scale validates its margins", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = product(foo(0 + member))),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) scale margin: foo"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = product(us(0 + member), us(0 + member))),
                data = dd, doFit = FALSE),
        "scale margins must be unique"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = product(ar1(0 + time))),
                data = dd, doFit = FALSE),
        "correlation-only"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = product(cs(0 + item))),
                data = dd, doFit = FALSE),
        "must match one of the specified margins"
    )
})

test_that("separable product syntax supports multi-column dense margins", {
    dd <- make_sep_dat()
    dd$x <- seq_len(nrow(dd))

    fit <- glmmTMB(y ~ 1 +
                       separable(us(0 + member + member:x) %x% ar1(0 + time) | group,
                                 scale = us(0 + member + member:x)),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepDims, c(4L, 3L))
    expect_equal(unname(fit$condReStruc[[1]]$blockSize), 12)
    expect_equal(fit$condReStruc[[1]]$blockNumTheta, 11)
    expect_equal(length(fit$condList$reTrms$cnms[[1]]), 12)
    expect_equal(head(fit$condList$reTrms$cnms[[1]], 4),
                 c("memberm1:time1", "memberm2:time1",
                   "memberm1:x:time1", "memberm2:x:time1"))
    expect_s3_class(fit$condList$reXterms[[1]], "separable_reXterms")
    expect_equal(fit$condList$reXterms[[1]]$margins$struc, c("us", "ar1"))
    expect_equal(length(fit$condList$reXterms[[1]]$terms), 2)
})

test_that("unsupported separable products fail at backend boundary", {
    dd <- expand.grid(member = factor(paste0("m", 1:2)),
                      time = factor(1:3),
                      item = factor(paste0("i", 1:2)),
                      group = factor(1:2))
    dd$y <- seq_len(nrow(dd)) / nrow(dd)

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(diag(0 + member) %x% ar1(0 + time) %x%
                                  homtoep(0 + item) | group,
                              scale = global()),
                data = dd, doFit = FALSE),
        "frontend parsed .*backend currently only evaluates"
    )
})

test_that("separable product margins must be no-intercept formulas", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(member) %x% ar1(0 + time) | group,
                              scale = us(member)),
                data = dd, doFit = FALSE),
        "no-intercept"
    )
})

test_that("separable product terms stay aligned after smooth augmentation", {
    skip_if_not_installed("mgcv")
    s <- mgcv::s
    dd <- make_sep_dat()
    dd$x <- seq_len(nrow(dd)) / nrow(dd)

    fit <- glmmTMB(y ~ s(x, k = 4) +
                       separable(us(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condList$ss, c("homdiag", "separable"))
    expect_equal(fit$condReStruc[[2]]$sepDims, c(2L, 3L))
    expect_s3_class(fit$condList$reXterms[[2]], "separable_reXterms")
})

test_that("separable supports explicit scale margin selection", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1 +
                       separable(ar1(0 + time) %x% us(0 + member) | group,
                                 scale = us(0 + member)),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[["us"]])))
    expect_equal(fit$condReStruc[[1]]$sepScaleSpec, 1L)
    expect_equal(fit$condReStruc[[1]]$blockNumTheta, 4)
})

test_that("multiple separable terms keep their spec order", {
    dd <- expand.grid(member = factor(c("A", "B")),
                      time = factor(1:3),
                      group1 = factor(1:2),
                      group2 = factor(1:2))
    dd$y <- 0

    fit <- glmmTMB(y ~ 1 +
                       (1 | group1) +
                       separable(us(0 + member) %x% ar1(0 + time) | group1) +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group2),
                   data = dd, doFit = FALSE)

    sep_terms <- vapply(fit$condReStruc, function(x) {
        identical(unname(x$blockCode), unname(.valid_covstruct[["separable"]]))
    }, logical(1))

    expect_equal(unname(which(sep_terms)), c(2L, 3L))
    expect_equal(fit$condReStruc[[2]]$blockNumTheta, 4)
    expect_equal(fit$condReStruc[[2]]$sepCodes,
                 unname(c(.valid_covstruct[["us"]], .valid_covstruct[["ar1"]])))
    expect_equal(fit$condReStruc[[3]]$blockNumTheta, 3)
    expect_equal(fit$condReStruc[[3]]$sepCodes,
                 unname(c(.valid_covstruct[["homcs"]], .valid_covstruct[["ar1"]])))
})

test_that("separable preserves the user-facing formula in the stored call", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    ftxt <- paste(deparse(fit$call$formula), collapse = " ")
    expect_match(ftxt, "homcs\\(0 \\+ member\\)")
    expect_false(grepl("data.frame", ftxt, fixed = TRUE))
})

test_that("separable preserves user-facing zi and dispersion formulas", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1,
                   ziformula = ~ separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   dispformula = ~ separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    ztxt <- paste(deparse(fit$call$ziformula), collapse = " ")
    dtxt <- paste(deparse(fit$call$dispformula), collapse = " ")
    expect_match(ztxt, "homcs\\(0 \\+ member\\)")
    expect_match(dtxt, "homcs\\(0 \\+ member\\)")
    expect_false(grepl("data.frame", ztxt, fixed = TRUE))
    expect_false(grepl("data.frame", dtxt, fixed = TRUE))
})

test_that("separable rejects unsupported margins", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% homcs(0 + time) | group),
                data = dd, doFit = FALSE),
        "specify the scale mode"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = ar1(0 + time)),
                data = dd, doFit = FALSE),
        "correlation-only"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(ar1(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "scale = global"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(ar1(0 + member) %x% ar1(0 + time) | group,
                              scale = global()),
                data = dd, doFit = FALSE),
        "backend currently only evaluates"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(ar1(0 + member) %x% ar1(0 + time) | group,
                              scale = product()),
                data = dd, doFit = FALSE),
        "needs at least one scale-capable margin"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(foo(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) margin: foo"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(un(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) margin: un"
    )
})

test_that("separable reports kronecker covariance for supported dense x ar1 pairs", {
    cases <- list(
        make_sep_case("cs", n_member = 5, n_time = 4),
        make_sep_case("homcs", n_member = 5, n_time = 4),
        make_sep_case("us", n_member = 4, n_time = 4),
        make_sep_case("us", reversed = TRUE),
        make_sep_case("cs", scale_mode = "global"),
        make_sep_case("homcs", scale_mode = "product"),
        make_sep_case("us", reversed = TRUE, scale_mode = "selected_product")
    )
    invisible(lapply(cases, expect_separable_case_vc))
})

test_that("separable likelihood matches dense MVN for supported dense x ar1 pairs", {
    cases <- list(
        make_sep_case("cs"),
        make_sep_case("cs", reversed = TRUE),
        make_sep_case("homcs"),
        make_sep_case("homcs", reversed = TRUE),
        make_sep_case("us", n_member = 3),
        make_sep_case("us", reversed = TRUE),
        make_sep_case("cs", scale_mode = "global"),
        make_sep_case("homcs", reversed = TRUE, scale_mode = "product"),
        make_sep_case("us", scale_mode = "selected_product")
    )
    invisible(lapply(cases, expect_separable_case_nll))
})

test_that("separable reports kronecker covariance for supported dense x dense pairs", {
    cases <- list(
        make_sep_dense_dense_case("cs", "cs", scale_mode = "global"),
        make_sep_dense_dense_case("homcs", "cs", scale_mode = "product"),
        make_sep_dense_dense_case("us", "homcs", scale_mode = "selected_first"),
        make_sep_dense_dense_case("cs", "us", scale_mode = "selected_second"),
        make_sep_dense_dense_case("us", "us", scale_mode = "product")
    )
    invisible(lapply(cases, expect_separable_case_vc))
})

test_that("separable likelihood matches dense MVN for supported dense x dense pairs", {
    cases <- list(
        make_sep_dense_dense_case("cs", "cs", scale_mode = "global"),
        make_sep_dense_dense_case("homcs", "cs", scale_mode = "product"),
        make_sep_dense_dense_case("us", "homcs", scale_mode = "selected_first"),
        make_sep_dense_dense_case("cs", "us", scale_mode = "selected_second"),
        make_sep_dense_dense_case("us", "us", scale_mode = "product")
    )
    invisible(lapply(cases, expect_separable_case_nll))
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
    R_member <- matrix(rho, n_member, n_member)
    diag(R_member) <- 1
    R_time <- outer(seq_len(n_time), seq_len(n_time),
                    function(i, j) phi^abs(i - j))
    R_full <- kronecker(R_time, R_member)
    sd_full <- rep(sd, n_time)
    Sigma <- diag(sd_full) %*% R_full %*% diag(sd_full)
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

test_that("separable simulation reports current limitation", {
    dd <- make_sep_dat(n_time = 2, reps = TRUE)

    theta <- c(log(1), qlogis((0.2 + 1) / 2), ar1_to_theta(0.3))
    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    expect_error(simulate(fit, nsim = 1),
                 "simulation is not yet implemented for separable covariance structures")
})
