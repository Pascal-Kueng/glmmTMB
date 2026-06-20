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
                                         "selected_product"),
                          ar1_struc = c("ar1", "hetar1")) {
    struc <- match.arg(struc)
    scale_mode <- match.arg(scale_mode)
    ar1_struc <- match.arg(ar1_struc)
    rho <- switch(struc, cs = 0.2, homcs = 0.3, us = -0.25)
    phi <- switch(struc, cs = 0.45, homcs = 0.4, us = 0.5)
    sd <- switch(struc,
        cs = seq(0.8, 1.2, length.out = n_member),
        homcs = 2,
        us = seq(0.8, 1.2, length.out = n_member)
    )
    time_sd <- if (ar1_struc == "hetar1") seq(1.1, 1.4, length.out = n_time)
               else rep(1.25, n_time)
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
    time_scale_theta <- if (scale_mode == "product") {
        if (ar1_struc == "hetar1") log(time_sd) else log(time_sd[[1]])
    } else numeric()
    time_theta <- c(time_scale_theta, ar1_to_theta(phi))
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
    ar1_call <- margin_call(ar1_struc, "time")
    member_sd_active <- if (scale_mode == "global") rep(1, n_member)
                        else if (length(sd) == 1) rep(sd, n_member) else sd
    time_sd_active <- if (scale_mode == "product") time_sd else rep(1, n_time)

    if (reversed) {
        form <- as.formula(as.call(list(quote(`~`), quote(y),
            as.call(list(quote(`+`), 1,
                         sep_call(as.call(list(as.name("%x%"),
                                               ar1_call, dense_call))))))))
        dense_form <- y ~ 1 + us(sepgrid(time, member) + 0 | group)
        theta <- if (scale_mode == "global") {
            c(log(global_sd), time_theta, member_theta)
        } else {
            c(time_theta, member_theta)
        }
        R_full <- kronecker(R_member, R_time)
        sd_full <- if (scale_mode == "global") {
            rep(global_sd, n_member * n_time)
        } else {
            as.vector(outer(time_sd_active, member_sd_active))
        }
        codes <- unname(c(.valid_covstruct[[ar1_struc]],
                          .valid_covstruct[[struc]]))
        kinds <- c(2L, 1L)
        scale_kinds <- c(if (ar1_struc == "hetar1") 2L else 1L,
                         if (struc == "homcs") 1L else 2L)
        scale_spec <- switch(scale_mode,
            global = integer(),
            margin = 1L,
            product = c(0L, 1L),
            selected_product = 1L)
    } else {
        form <- as.formula(as.call(list(quote(`~`), quote(y),
            as.call(list(quote(`+`), 1,
                         sep_call(as.call(list(as.name("%x%"),
                                               dense_call, ar1_call))))))))
        dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
        theta <- if (scale_mode == "global") {
            c(log(global_sd), member_theta, time_theta)
        } else {
            c(member_theta, time_theta)
        }
        R_full <- kronecker(R_time, R_member)
        sd_full <- if (scale_mode == "global") {
            rep(global_sd, n_member * n_time)
        } else {
            as.vector(outer(member_sd_active, time_sd_active))
        }
        codes <- unname(c(.valid_covstruct[[struc]],
                          .valid_covstruct[[ar1_struc]]))
        kinds <- c(1L, 2L)
        scale_kinds <- c(if (struc == "homcs") 1L else 2L,
                         if (ar1_struc == "hetar1") 2L else 1L)
        scale_spec <- switch(scale_mode,
            global = integer(),
            margin = 0L,
            product = c(0L, 1L),
            selected_product = 0L)
    }

    scale_mode_code <- switch(scale_mode,
        margin = 1L, global = 2L, product = 3L, selected_product = 4L)
    list(form = form, dense_form = dense_form, theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full, sd_full = sd_full,
         codes = codes, kinds = kinds, dispatch = 1L,
         scale_kinds = scale_kinds,
         scale_mode = scale_mode_code,
         scale_spec = scale_spec,
         n_member = n_member, n_time = n_time)
}

case_data <- function(case, reps = FALSE, n_group = 2) {
    if (!is.null(case$dd)) return(case$dd)
    make_sep_dat(n_member = case$n_member, n_time = case$n_time,
                 n_group = n_group, reps = reps)
}

expected_sep_scale_kinds <- function(case) {
    if (!is.null(case$scale_kinds)) return(case$scale_kinds)
    struc <- names(.valid_covstruct)[match(case$codes, unname(.valid_covstruct))]
    ifelse(struc %in% c("homdiag", "homcs", "ar1", "ou", "exp",
                        "gau", "mat", "homtoep"), 1L, 2L)
}

expect_separable_case_vc <- function(case) {
    dd <- case_data(case, reps = TRUE)
    fit <- fit_fixed_theta(case$form, dd, case$theta)
    restruc <- fit$modelInfo$reStruc$condReStruc[[1]]
    vc <- VarCorr(fit)$cond[[1]]

    expect_equal(restruc$sepCodes, case$codes)
    expect_equal(restruc$sepBuilderKinds, case$kinds)
    expect_equal(restruc$sepScaleKinds, expected_sep_scale_kinds(case))
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

make_diag_margin <- function(struc = c("diag", "homdiag"), n = 2,
                             sd = seq(0.7, 1.1, length.out = n)) {
    struc <- match.arg(struc)
    list(
        struc = struc,
        sd = if (struc == "homdiag") sd[[1]] else sd,
        R = diag(n),
        scale_theta = log(if (struc == "homdiag") sd[[1]] else sd),
        corr_theta = numeric()
    )
}

make_ar1_margin <- function(struc = c("ar1", "hetar1"), n = 3,
                            sd = seq(0.9, 1.3, length.out = n),
                            phi = 0.4) {
    struc <- match.arg(struc)
    R <- outer(seq_len(n), seq_len(n), function(i, j) phi^abs(i - j))
    list(
        struc = struc,
        sd = if (struc == "hetar1") sd else rep(sd[[1]], n),
        R = R,
        scale_theta = if (struc == "hetar1") log(sd) else log(sd[[1]]),
        corr_theta = ar1_to_theta(phi)
    )
}

make_toep_margin <- function(struc = c("toep", "homtoep"), n = 3,
                             sd = seq(0.8, 1.2, length.out = n),
                             rho = seq(0.25, 0.05, length.out = n - 1L)) {
    struc <- match.arg(struc)
    R <- diag(n)
    for (lag in seq_len(n - 1L)) {
        R[row(R) == col(R) + lag] <- rho[[lag]]
        R[col(R) == row(R) + lag] <- rho[[lag]]
    }
    list(
        struc = struc,
        sd = if (struc == "homtoep") sd[[1]] else sd,
        R = R,
        scale_theta = log(if (struc == "homtoep") sd[[1]] else sd),
        corr_theta = ar1_to_theta(rho)
    )
}

make_spatial_margin <- function(struc = c("ou", "exp", "gau", "mat"),
                                coords = cbind(seq_len(3)), theta = 0.2,
                                sd = 1.15) {
    struc <- match.arg(struc)
    coords <- as.matrix(coords)
    D <- as.matrix(dist(coords))
    R <- switch(struc,
        ou = exp(-exp(theta) * D),
        exp = exp(-D * exp(-theta)),
        gau = exp(-(D^2) * exp(-2 * theta)),
        mat = {
            phi <- exp(theta[[1]])
            kappa <- exp(theta[[2]])
            M <- matrix(1, nrow(D), ncol(D))
            keep <- D > 0
            x <- D[keep] / phi
            M[keep] <- x^kappa * besselK(x, kappa) /
                (gamma(kappa) * 2^(kappa - 1))
            M
        }
    )
    diag(R) <- 1
    list(
        struc = struc,
        sd = sd,
        R = R,
        scale_theta = log(sd),
        corr_theta = theta,
        levels = paste0("(", apply(coords, 1, paste, collapse = ","), ")")
    )
}

margin_call <- function(struc, var) {
    as.call(list(as.name(struc),
                 as.call(list(as.name("+"), 0, as.name(var)))))
}

sep_kind_code <- function(struc) {
    switch(struc,
           cs = 1L, homcs = 1L, us = 1L,
           ar1 = 2L, hetar1 = 2L,
           diag = 3L, homdiag = 3L,
           ou = 4L, exp = 4L, gau = 4L, mat = 4L,
           toep = 5L, homtoep = 5L,
           propto = 6L, equalto = 6L)
}

make_sep_margin_pair_case <- function(m0, m1,
                                      var0 = "member", var1 = "item",
                                      scale_mode = c("global", "product",
                                                     "selected_first",
                                                     "selected_second")) {
    scale_mode <- match.arg(scale_mode)
    n0 <- nrow(m0$R)
    n1 <- nrow(m1$R)
    global_sd <- 1.25
    dd <- expand.grid(factor(seq_len(n0)), factor(seq_len(n1)),
                      group = factor(seq_len(2)))
    names(dd)[1:2] <- c(var0, var1)
    dd$y <- 0

    call0 <- margin_call(m0$struc, var0)
    call1 <- margin_call(m1$struc, var1)
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

    scale_index <- switch(scale_mode,
        global = integer(),
        product = which(c(length(m0$scale_theta), length(m1$scale_theta)) > 0L),
        selected_first = 1L,
        selected_second = 2L
    )
    theta0 <- c(if (1L %in% scale_index) m0$scale_theta, m0$corr_theta)
    theta1 <- c(if (2L %in% scale_index) m1$scale_theta, m1$corr_theta)
    theta <- if (scale_mode == "global") c(log(global_sd), theta0, theta1)
             else c(theta0, theta1)

    sd0 <- if (1L %in% scale_index) {
        if (length(m0$sd) == 1L) rep(m0$sd, n0) else m0$sd
    } else rep(1, n0)
    sd1 <- if (2L %in% scale_index) {
        if (length(m1$sd) == 1L) rep(m1$sd, n1) else m1$sd
    } else rep(1, n1)
    sd_full <- if (scale_mode == "global") rep(global_sd, n0 * n1)
               else as.vector(outer(sd0, sd1))
    R_full <- kronecker(m1$R, m0$R)

    dense_rhs <- as.call(list(
        as.name("|"),
        as.call(list(as.name("+"),
                     as.call(list(as.name("sepgrid"),
                                  as.name(var0), as.name(var1))),
                     0)),
        quote(group)
    ))
    dense_form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1,
                     as.call(list(as.name("us"), dense_rhs)))))))

    list(form = form,
         dense_form = dense_form,
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = unname(c(.valid_covstruct[[m0$struc]],
                          .valid_covstruct[[m1$struc]])),
         kinds = c(sep_kind_code(m0$struc), sep_kind_code(m1$struc)),
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             global = 2L, product = 3L,
             selected_first = 4L, selected_second = 4L),
         scale_spec = as.integer(scale_index - 1L))
}

make_sep_margin_chain_case <- function(margins, vars = paste0("v", seq_along(margins)),
                                       scale_mode = c("global", "product",
                                                      "selected_product"),
                                       scale_index = NULL) {
    scale_mode <- match.arg(scale_mode)
    dims <- vapply(margins, function(m) nrow(m$R), integer(1))
    factors <- Map(function(m, n) {
        if (!is.null(m$levels)) factor(m$levels, levels = m$levels)
        else factor(seq_len(n))
    }, margins, dims)
    names(factors) <- vars
    dd <- do.call(expand.grid, c(factors, list(group = factor(seq_len(2)))))
    dd$y <- 0

    calls <- Map(margin_call, vapply(margins, `[[`, character(1), "struc"),
                 vars)
    product_call <- Reduce(function(a, b) {
        as.call(list(as.name("%x%"), a, b))
    }, calls)
    if (scale_mode == "global") {
        scale_call <- quote(global())
        scale_index <- integer()
    } else if (scale_mode == "product") {
        scale_call <- quote(product())
        scale_index <- which(vapply(margins, function(m) length(m$scale_theta),
                                    integer(1)) > 0L)
    } else {
        scale_call <- as.call(c(list(as.name("product")), calls[scale_index]))
    }
    sep_call <- as.call(list(
        as.name("separable"),
        as.call(list(as.name("|"), product_call, quote(group))),
        scale = scale_call
    ))
    form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1, sep_call)))))

    theta <- unlist(lapply(seq_along(margins), function(i) {
        c(if (i %in% scale_index) margins[[i]]$scale_theta,
          margins[[i]]$corr_theta)
    }), use.names = FALSE)
    global_sd <- 1.25
    if (scale_mode == "global") theta <- c(log(global_sd), theta)

    coords <- do.call(expand.grid, lapply(dims, seq_len))
    sd_full <- if (scale_mode == "global") {
        rep(global_sd, prod(dims))
    } else {
        ans <- rep(1, prod(dims))
        for (i in scale_index) {
            sd_i <- margins[[i]]$sd
            if (length(sd_i) == 1L) sd_i <- rep(sd_i, dims[[i]])
            ans <- ans * sd_i[coords[[i]]]
        }
        ans
    }
    R_full <- margins[[1]]$R
    for (i in seq_along(margins)[-1]) {
        R_full <- kronecker(margins[[i]]$R, R_full)
    }

    dense_grid <- as.call(c(list(as.name("sepgrid")), lapply(vars, as.name)))
    dense_rhs <- as.call(list(
        as.name("|"),
        as.call(list(as.name("+"), dense_grid, 0)),
        quote(group)
    ))
    dense_form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1,
                     as.call(list(as.name("us"), dense_rhs)))))))

    list(form = form,
         dense_form = dense_form,
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = unname(vapply(margins, function(m) .valid_covstruct[[m$struc]],
                               numeric(1))),
         kinds = vapply(margins, function(m) sep_kind_code(m$struc),
                        integer(1)),
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             global = 2L, product = 3L, selected_product = 4L),
         scale_spec = as.integer(scale_index - 1L))
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
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             global = 2L, product = 3L,
             selected_first = 4L, selected_second = 4L),
         scale_spec = switch(scale_mode,
             global = integer(), product = c(0L, 1L),
             selected_first = 0L, selected_second = 1L))
}

make_sep_ar1_ar1_case <- function(struc0 = c("ar1", "hetar1"),
                                  struc1 = c("ar1", "hetar1"),
                                  scale_mode = c("global", "margin",
                                                 "product",
                                                 "selected_first",
                                                 "selected_second"),
                                  n0 = 3, n1 = 4, phi0 = 0.35, phi1 = 0.55) {
    struc0 <- match.arg(struc0)
    struc1 <- match.arg(struc1)
    scale_mode <- match.arg(scale_mode)
    m0 <- make_ar1_margin(struc0, n0, phi = phi0)
    m1 <- make_ar1_margin(struc1, n1,
                          sd = seq(1.1, 1.5, length.out = n1),
                          phi = phi1)
    global_sd <- 1.25
    dd <- expand.grid(time0 = factor(seq_len(n0)),
                      time1 = factor(seq_len(n1)),
                      group = factor(seq_len(2)))
    dd$y <- 0

    call0 <- margin_call(struc0, "time0")
    call1 <- margin_call(struc1, "time1")
    scale_call <- switch(scale_mode,
        global = quote(global()),
        margin = NULL,
        product = quote(product()),
        selected_first = as.call(list(as.name("product"), call0)),
        selected_second = as.call(list(as.name("product"), call1))
    )
    sep_args <- list(
        as.name("separable"),
        as.call(list(as.name("|"),
                     as.call(list(as.name("%x%"), call0, call1)),
                     quote(group)))
    )
    if (!is.null(scale_call)) sep_args$scale <- scale_call
    sep_call <- as.call(sep_args)
    form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1, sep_call)))))

    scale_index <- switch(scale_mode,
        global = integer(),
        margin = which(c(struc0, struc1) == "hetar1"),
        product = seq_along(c(struc0, struc1)),
        selected_first = 1L,
        selected_second = 2L)
    theta0 <- c(if (1L %in% scale_index) m0$scale_theta, m0$corr_theta)
    theta1 <- c(if (2L %in% scale_index) m1$scale_theta, m1$corr_theta)
    theta <- if (scale_mode == "global") c(log(global_sd), theta0, theta1)
             else c(theta0, theta1)

    R_full <- kronecker(m1$R, m0$R)
    sd0 <- if (1L %in% scale_index) m0$sd else rep(1, n0)
    sd1 <- if (2L %in% scale_index) m1$sd else rep(1, n1)
    sd_full <- if (scale_mode == "global") rep(global_sd, n0 * n1)
               else as.vector(outer(sd0, sd1))

    list(form = form,
         dense_form = y ~ 1 + us(sepgrid(time0, time1) + 0 | group),
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = unname(c(.valid_covstruct[[struc0]],
                          .valid_covstruct[[struc1]])),
         kinds = c(2L, 2L),
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             margin = 1L, global = 2L, product = 3L,
             selected_first = 4L, selected_second = 4L),
         scale_spec = as.integer(scale_index - 1L))
}

make_sep_diag_ar1_case <- function(struc = c("diag", "homdiag"),
                                   reversed = FALSE,
                                   scale_mode = c("margin", "global",
                                                  "product", "selected"),
                                   n_diag = 3, n_time = 4, phi = 0.45,
                                   ar1_struc = c("ar1", "hetar1")) {
    struc <- match.arg(struc)
    scale_mode <- match.arg(scale_mode)
    ar1_struc <- match.arg(ar1_struc)
    m <- make_diag_margin(struc, n_diag)
    a <- make_ar1_margin(ar1_struc, n_time, phi = phi)
    global_sd <- 1.2
    dd <- expand.grid(member = factor(paste0("m", seq_len(n_diag))),
                      time = factor(seq_len(n_time)),
                      group = factor(seq_len(2)))
    dd$y <- 0

    diag_call <- margin_call(struc, "member")
    ar1_call <- margin_call(ar1_struc, "time")
    scale_call <- switch(scale_mode,
        margin = NULL,
        global = quote(global()),
        product = quote(product()),
        selected = as.call(list(as.name("product"), diag_call))
    )
    lhs <- if (reversed) as.call(list(as.name("%x%"), ar1_call, diag_call))
           else as.call(list(as.name("%x%"), diag_call, ar1_call))
    sep_call <- as.call(c(list(as.name("separable"),
                               as.call(list(as.name("|"), lhs, quote(group)))),
                          if (is.null(scale_call)) list()
                          else list(scale = scale_call)))
    form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1, sep_call)))))

    sd_diag <- if (length(m$sd) == 1L) rep(m$sd, n_diag) else m$sd
    sd_time <- a$sd
    diag_active <- scale_mode != "global"
    time_active <- scale_mode == "product"
    theta_diag <- if (diag_active) m$scale_theta else numeric()
    theta_time <- c(if (time_active) a$scale_theta, a$corr_theta)
    theta <- if (reversed) {
        c(if (scale_mode == "global") log(global_sd),
          theta_time, theta_diag)
    } else {
        c(if (scale_mode == "global") log(global_sd),
          theta_diag, theta_time)
    }
    sd_full <- if (scale_mode == "global") {
        rep(global_sd, n_diag * n_time)
    } else if (reversed) {
        as.vector(outer(if (time_active) sd_time else rep(1, n_time),
                        sd_diag))
    } else {
        as.vector(outer(sd_diag,
                        if (time_active) sd_time else rep(1, n_time)))
    }
    R_full <- if (reversed) kronecker(m$R, a$R)
              else kronecker(a$R, m$R)
    scale_spec <- switch(scale_mode,
        global = integer(),
        margin = if (reversed) 1L else 0L,
        product = c(0L, 1L),
        selected = if (reversed) 1L else 0L)

    list(form = form,
         dense_form = if (reversed) y ~ 1 + us(sepgrid(time, member) + 0 | group)
                      else y ~ 1 + us(sepgrid(member, time) + 0 | group),
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = if (reversed) unname(c(.valid_covstruct[[ar1_struc]],
                                        .valid_covstruct[[struc]]))
                 else unname(c(.valid_covstruct[[struc]],
                               .valid_covstruct[[ar1_struc]])),
         kinds = if (reversed) c(2L, 3L) else c(3L, 2L),
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             margin = 1L, global = 2L, product = 3L, selected = 4L),
         scale_spec = scale_spec)
}

make_sep_diag_dense_case <- function(diag_struc = c("diag", "homdiag"),
                                     dense_struc = c("cs", "homcs", "us"),
                                     reversed = FALSE,
                                     scale_mode = c("global", "product",
                                                    "selected_diag",
                                                    "selected_dense"),
                                     n_diag = 3, n_dense = 3) {
    diag_struc <- match.arg(diag_struc)
    dense_struc <- match.arg(dense_struc)
    scale_mode <- match.arg(scale_mode)
    d <- make_diag_margin(diag_struc, n_diag)
    m <- make_dense_margin(dense_struc, n_dense, rho = 0.25)
    global_sd <- 1.15
    dd <- expand.grid(member = factor(paste0("m", seq_len(n_diag))),
                      item = factor(paste0("i", seq_len(n_dense))),
                      group = factor(seq_len(2)))
    dd$y <- 0

    diag_call <- margin_call(diag_struc, "member")
    dense_call <- margin_call(dense_struc, "item")
    scale_call <- switch(scale_mode,
        global = quote(global()),
        product = quote(product()),
        selected_diag = as.call(list(as.name("product"), diag_call)),
        selected_dense = as.call(list(as.name("product"), dense_call))
    )
    lhs <- if (reversed) as.call(list(as.name("%x%"), dense_call, diag_call))
           else as.call(list(as.name("%x%"), diag_call, dense_call))
    sep_call <- as.call(list(
        as.name("separable"),
        as.call(list(as.name("|"), lhs, quote(group))),
        scale = scale_call
    ))
    form <- as.formula(as.call(list(quote(`~`), quote(y),
        as.call(list(quote(`+`), 1, sep_call)))))

    theta_diag <- switch(scale_mode,
        global = numeric(),
        product = d$scale_theta,
        selected_diag = d$scale_theta,
        selected_dense = numeric()
    )
    theta_dense <- switch(scale_mode,
        global = m$corr_theta,
        product = c(m$scale_theta, m$corr_theta),
        selected_diag = m$corr_theta,
        selected_dense = c(m$scale_theta, m$corr_theta)
    )
    theta <- if (reversed) {
        c(if (scale_mode == "global") log(global_sd), theta_dense, theta_diag)
    } else {
        c(if (scale_mode == "global") log(global_sd), theta_diag, theta_dense)
    }

    sd_diag <- if (length(d$sd) == 1L) rep(d$sd, n_diag) else d$sd
    sd_dense <- if (length(m$sd) == 1L) rep(m$sd, n_dense) else m$sd
    sd_full <- switch(scale_mode,
        global = rep(global_sd, n_diag * n_dense),
        product = if (reversed) as.vector(outer(sd_dense, sd_diag))
                  else as.vector(outer(sd_diag, sd_dense)),
        selected_diag = if (reversed) rep(sd_diag, each = n_dense)
                        else rep(sd_diag, n_dense),
        selected_dense = if (reversed) rep(sd_dense, n_diag)
                         else rep(sd_dense, each = n_diag)
    )
    R_full <- if (reversed) kronecker(d$R, m$R) else kronecker(m$R, d$R)

    list(form = form,
         dense_form = if (reversed) y ~ 1 + us(sepgrid(item, member) + 0 | group)
                      else y ~ 1 + us(sepgrid(member, item) + 0 | group),
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = if (reversed) unname(c(.valid_covstruct[[dense_struc]],
                                        .valid_covstruct[[diag_struc]]))
                 else unname(c(.valid_covstruct[[diag_struc]],
                               .valid_covstruct[[dense_struc]])),
         kinds = if (reversed) c(1L, 3L) else c(3L, 1L),
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             global = 2L, product = 3L,
             selected_diag = 4L, selected_dense = 4L),
         scale_spec = switch(scale_mode,
             global = integer(),
             product = c(0L, 1L),
             selected_diag = if (reversed) 1L else 0L,
             selected_dense = if (reversed) 0L else 1L))
}

make_sep_diag_diag_case <- function(struc0 = c("diag", "homdiag"),
                                    struc1 = c("diag", "homdiag"),
                                    scale_mode = c("global", "product",
                                                   "selected_first"),
                                    n0 = 2, n1 = 3) {
    struc0 <- match.arg(struc0)
    struc1 <- match.arg(struc1)
    scale_mode <- match.arg(scale_mode)
    m0 <- make_diag_margin(struc0, n0)
    m1 <- make_diag_margin(struc1, n1, sd = seq(1.1, 1.4, length.out = n1))
    global_sd <- 1.3
    dd <- expand.grid(member = factor(paste0("m", seq_len(n0))),
                      item = factor(paste0("i", seq_len(n1))),
                      group = factor(seq_len(2)))
    dd$y <- 0

    call0 <- margin_call(struc0, "member")
    call1 <- margin_call(struc1, "item")
    scale_call <- switch(scale_mode,
        global = quote(global()),
        product = quote(product()),
        selected_first = as.call(list(as.name("product"), call0))
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
        global = numeric(),
        product = m0$scale_theta,
        selected_first = m0$scale_theta
    )
    theta1 <- switch(scale_mode,
        global = numeric(),
        product = m1$scale_theta,
        selected_first = numeric()
    )
    theta <- c(if (scale_mode == "global") log(global_sd), theta0, theta1)
    sd0 <- if (length(m0$sd) == 1L) rep(m0$sd, n0) else m0$sd
    sd1 <- if (length(m1$sd) == 1L) rep(m1$sd, n1) else m1$sd
    sd_full <- switch(scale_mode,
        global = rep(global_sd, n0 * n1),
        product = as.vector(outer(sd0, sd1)),
        selected_first = rep(sd0, n1)
    )
    R_full <- diag(n0 * n1)

    list(form = form,
         dense_form = y ~ 1 + us(sepgrid(member, item) + 0 | group),
         dd = dd,
         theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full,
         sd_full = sd_full,
         codes = unname(c(.valid_covstruct[[struc0]],
                          .valid_covstruct[[struc1]])),
         kinds = c(3L, 3L),
         dispatch = 1L,
         scale_mode = switch(scale_mode,
             global = 2L, product = 3L, selected_first = 4L),
         scale_spec = switch(scale_mode,
             global = integer(), product = c(0L, 1L), selected_first = 0L))
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
    expect_equal(h$condReStruc[[1]]$sepBuilderKinds, c(2L, 1L))
    expect_equal(h$condReStruc[[1]]$sepScaleKinds, c(1L, 1L))
    expect_equal(h$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleSpec, 1L)

    expect_equal(u$condReStruc[[1]]$blockNumTheta, 4)
    expect_equal(u$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["us"]], .valid_covstruct[["ar1"]])))
    expect_equal(u$condReStruc[[1]]$sepBuilderKinds, c(1L, 2L))
    expect_equal(u$condReStruc[[1]]$sepScaleKinds, c(2L, 1L))
    expect_equal(u$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(u$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(u$condReStruc[[1]]$sepScaleSpec, 0L)
})

test_that("separable parser builds structured specs from splitForm output", {
    f <- y ~ 1 +
        separable(homcs(0 + member) %x% ar1(0 + time) | group,
                  scale = homcs(0 + member))
    ss <- reformulas::splitForm(f, specials = c(names(.valid_covstruct), "s"))
    specs <- glmmTMB:::.sep_specs_from_split(ss)
    spec <- specs[[1]]

    expect_length(specs, 1)
    expect_equal(spec$grid, c("member", "time"))
    expect_equal(unname(spec$margins$struc), c("homcs", "ar1"))
    expect_equal(spec$scale$mode, "margin")
    expect_equal(unname(spec$scale$margins$struc), "homcs")
    expect_equal(unname(spec$scale$margins$var), "member")
    expect_equal(ss$reTrmClasses, "separable")
    expect_equal(deparse(ss$reTrmFormulas[[1]]),
                 "homcs(0 + member) %x% ar1(0 + time) | group")
    expect_equal(deparse(ss$reTrmAddArgs[[1]]),
                 "separable(scale = homcs(0 + member))")
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

test_that("separable frontend parses simple existing covariance margins", {
    f <- y ~ 1 +
        separable(diag(0 + member) %x% ar1(0 + time) %x%
                      homtoep(0 + item) | group,
                  scale = product(diag(0 + member), homtoep(0 + item)))
    ss <- reformulas::splitForm(f, specials = c(names(.valid_covstruct), "s"))
    spec <- glmmTMB:::.sep_specs_from_split(ss)[[1]]

    expect_equal(spec$grid, c("member", "time", "item"))
    expect_equal(spec$margins$struc, c("diag", "ar1", "homtoep"))
    expect_equal(spec$scale$mode, "selected_product")
    expect_equal(spec$scale$margins$struc, c("diag", "homtoep"))
})

test_that("separable parser carries margin extra arguments", {
    m <- glmmTMB:::.sep_product_margin_spec(quote(ar1(0 + time, rho_source)))
    m2 <- glmmTMB:::.sep_product_margin_spec(quote(ar1(0 + time, other_source)))

    expect_equal(m$struc, "ar1")
    expect_equal(m$var, "time")
    expect_equal(length(m$extra[[1]]), 1L)
    expect_equal(deparse(m$extra[[1]][[1]]), "rho_source")
    expect_false(identical(glmmTMB:::.sep_margin_key(m),
                           glmmTMB:::.sep_margin_key(m2)))
})

test_that("separable parser records global and product scale modes", {
    f_global <- y ~ 1 +
        separable(ar1(0 + member) %x% ar1(0 + time) | group,
                  scale = global())
    ss_global <- reformulas::splitForm(f_global,
                                       specials = c(names(.valid_covstruct), "s"))
    expect_equal(glmmTMB:::.sep_specs_from_split(ss_global)[[1]]$scale$mode,
                 "global")

    f_product <- y ~ 1 +
        separable(us(0 + member) %x% cs(0 + item) | group,
                  scale = product())
    ss_product <- reformulas::splitForm(f_product,
                                        specials = c(names(.valid_covstruct), "s"))
    expect_equal(glmmTMB:::.sep_specs_from_split(ss_product)[[1]]$scale$mode,
                 "product")
})

test_that("separable records registry-derived theta blocks", {
    dd <- make_sep_dat()
    fit <- glmmTMB(y ~ 1 +
                       separable(ar1(0 + member) %x% ar1(0 + time) | group,
                                 scale = global()),
                   data = dd, doFit = FALSE)
    restruc <- fit$condReStruc[[1]]

    expect_equal(restruc$sepThetaBlockMargins, c(-1L, 0L, 1L))
    expect_equal(restruc$sepThetaBlockKinds, c(1L, 3L, 3L))
    expect_equal(restruc$sepThetaBlockStarts, 0:2)
    expect_equal(restruc$sepThetaBlockLengths, c(1L, 1L, 1L))
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
                              scale = product(ar1(0 + time, foo))),
                data = dd, doFit = FALSE),
        "scale margin ar1\\(time\\) takes 0 extra arguments, but got 1"
    )
    fit <- glmmTMB(y ~ 1 +
                       separable(us(0 + member) %x% ar1(0 + time) | group,
                                 scale = product(ar1(0 + time))),
                   data = dd, doFit = FALSE)
    expect_equal(fit$condReStruc[[1]]$sepScaleSpec, 1L)
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

test_that("separable spatial margins require numeric coordinate levels", {
    dd <- expand.grid(member = factor(paste0("m", 1:2)),
                      time = factor(1:3),
                      item = factor(paste0("i", 1:2)),
                      group = factor(1:2))
    dd$y <- seq_len(nrow(dd)) / nrow(dd)

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(diag(0 + member) %x% ou(0 + time) %x%
                                  homtoep(0 + item) | group,
                              scale = global()),
                data = dd, doFit = FALSE),
        "spatial margins require numeric coordinate levels"
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

test_that("separable validates unsupported margins and scale choices", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% homcs(0 + time) | group),
                data = dd, doFit = FALSE),
        "specify the scale mode"
    )
    fit <- glmmTMB(y ~ 1 +
                       separable(us(0 + member) %x% ar1(0 + time) | group,
                                 scale = ar1(0 + time)),
                   data = dd, doFit = FALSE)
    expect_equal(fit$condReStruc[[1]]$sepScaleSpec, 1L)
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(ar1(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "no unambiguous scale margin"
    )
    fit <- glmmTMB(y ~ 1 +
                       separable(ar1(0 + member) %x% ar1(0 + time) | group,
                                 scale = product()),
                   data = dd, doFit = FALSE)
    expect_equal(fit$condReStruc[[1]]$sepScaleSpec, c(0L, 1L))
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(foo(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) margin: foo"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member, extra) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "margin us\\(member\\) takes 0 extra arguments, but got 1"
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
        make_sep_case("us", reversed = TRUE, scale_mode = "selected_product"),
        make_sep_case("cs", scale_mode = "product", ar1_struc = "hetar1")
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
        make_sep_case("us", scale_mode = "selected_product"),
        make_sep_case("us", reversed = TRUE, scale_mode = "product",
                      ar1_struc = "hetar1")
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

test_that("separable reports kronecker covariance for ar1 x ar1", {
    cases <- list(
        make_sep_ar1_ar1_case(),
        make_sep_ar1_ar1_case("ar1", "ar1", scale_mode = "product"),
        make_sep_ar1_ar1_case("hetar1", "ar1", scale_mode = "margin"),
        make_sep_ar1_ar1_case("ar1", "hetar1",
                              scale_mode = "selected_second"),
        make_sep_ar1_ar1_case("hetar1", "hetar1", scale_mode = "product")
    )
    invisible(lapply(cases, expect_separable_case_vc))
})

test_that("separable likelihood matches dense MVN for ar1 x ar1", {
    cases <- list(
        make_sep_ar1_ar1_case(),
        make_sep_ar1_ar1_case("ar1", "ar1", scale_mode = "product"),
        make_sep_ar1_ar1_case("hetar1", "ar1", scale_mode = "margin"),
        make_sep_ar1_ar1_case("ar1", "hetar1",
                              scale_mode = "selected_second"),
        make_sep_ar1_ar1_case("hetar1", "hetar1", scale_mode = "product")
    )
    invisible(lapply(cases, expect_separable_case_nll))
})

test_that("separable reports kronecker covariance for diagonal margin pairs", {
    cases <- list(
        make_sep_diag_ar1_case("diag"),
        make_sep_diag_ar1_case("homdiag", reversed = TRUE,
                               scale_mode = "product"),
        make_sep_diag_ar1_case("diag", scale_mode = "product",
                               ar1_struc = "hetar1"),
        make_sep_diag_dense_case("diag", "us", scale_mode = "product"),
        make_sep_diag_dense_case("homdiag", "cs", reversed = TRUE,
                                 scale_mode = "selected_diag"),
        make_sep_diag_diag_case("diag", "homdiag", scale_mode = "product"),
        make_sep_diag_diag_case("homdiag", "diag", scale_mode = "global")
    )
    invisible(lapply(cases, expect_separable_case_vc))
})

test_that("separable likelihood matches dense MVN for diagonal margin pairs", {
    cases <- list(
        make_sep_diag_ar1_case("diag"),
        make_sep_diag_ar1_case("homdiag", reversed = TRUE,
                               scale_mode = "selected"),
        make_sep_diag_ar1_case("homdiag", reversed = TRUE,
                               scale_mode = "product",
                               ar1_struc = "hetar1"),
        make_sep_diag_dense_case("diag", "us", scale_mode = "product"),
        make_sep_diag_dense_case("homdiag", "homcs",
                                 scale_mode = "selected_dense"),
        make_sep_diag_diag_case("diag", "homdiag", scale_mode = "product"),
        make_sep_diag_diag_case("homdiag", "diag",
                                scale_mode = "selected_first")
    )
    invisible(lapply(cases, expect_separable_case_nll))
})

test_that("separable supports Toeplitz correlation margins", {
    cases <- list(
        make_sep_margin_pair_case(
            make_toep_margin("toep", n = 3),
            make_dense_margin("us", n = 3),
            scale_mode = "product"
        ),
        make_sep_margin_pair_case(
            make_toep_margin("homtoep", n = 4),
            make_ar1_margin("ar1", n = 3),
            var0 = "time0", var1 = "time1",
            scale_mode = "selected_first"
        ),
        make_sep_margin_pair_case(
            make_diag_margin("diag", n = 3),
            make_toep_margin("homtoep", n = 3),
            scale_mode = "global"
        )
    )
    invisible(lapply(cases, expect_separable_case_vc))
    invisible(lapply(cases, expect_separable_case_nll))
})

test_that("separable supports three to five correlation-matrix margins", {
    cases <- list(
        make_sep_margin_chain_case(
            list(make_dense_margin("us", n = 2),
                 make_ar1_margin("hetar1", n = 3),
                 make_toep_margin("homtoep", n = 2)),
            vars = c("member", "time", "item"),
            scale_mode = "product"
        ),
        make_sep_margin_chain_case(
            list(make_diag_margin("diag", n = 2),
                 make_ar1_margin("ar1", n = 3),
                 make_dense_margin("cs", n = 2)),
            vars = c("member", "time", "item"),
            scale_mode = "global"
        ),
        make_sep_margin_chain_case(
            list(make_dense_margin("homcs", n = 2),
                 make_diag_margin("homdiag", n = 3),
                 make_toep_margin("toep", n = 2)),
            vars = c("member", "item", "time"),
            scale_mode = "selected_product",
            scale_index = c(1L, 3L)
        ),
        make_sep_margin_chain_case(
            list(make_dense_margin("cs", n = 2),
                 make_diag_margin("diag", n = 2),
                 make_ar1_margin("ar1", n = 2),
                 make_toep_margin("homtoep", n = 2)),
            vars = c("member", "item", "time", "occasion"),
            scale_mode = "selected_product",
            scale_index = 2L
        ),
        make_sep_margin_chain_case(
            list(make_dense_margin("homcs", n = 2),
                 make_diag_margin("homdiag", n = 2),
                 make_ar1_margin("ar1", n = 2),
                 make_toep_margin("homtoep", n = 2),
                 make_dense_margin("cs", n = 2)),
            vars = c("member", "item", "time", "occasion", "rater"),
            scale_mode = "global"
        )
    )
    invisible(lapply(cases, expect_separable_case_vc))
    invisible(lapply(cases, expect_separable_case_nll))
})

test_that("separable supports spatial correlation margins", {
    cases <- list(
        make_sep_margin_chain_case(
            list(make_spatial_margin("ou", coords = cbind(c(0, 0.5, 2)),
                                     theta = -0.1),
                 make_dense_margin("us", n = 2)),
            vars = c("time", "member"),
            scale_mode = "global"
        ),
        make_sep_margin_chain_case(
            list(make_spatial_margin("exp", coords = cbind(c(0, 1, 2),
                                                           c(0, 0, 1)),
                                     theta = 0.3),
                 make_diag_margin("diag", n = 2),
                 make_spatial_margin("gau", coords = cbind(c(0, 2)),
                                     theta = 0.7)),
            vars = c("space", "member", "time"),
            scale_mode = "product"
        ),
        make_sep_margin_chain_case(
            list(make_spatial_margin("mat", coords = cbind(c(0, 1, 3)),
                                     theta = c(log(1.2), log(0.8))),
                 make_diag_margin("homdiag", n = 2)),
            vars = c("space", "member"),
            scale_mode = "global"
        )
    )
    invisible(lapply(cases, expect_separable_case_vc))
    invisible(lapply(cases, expect_separable_case_nll))
})

test_that("separable supports propto matrix margins", {
    n_member <- 3
    n_time <- 4
    dd <- expand.grid(member = factor(paste0("m", seq_len(n_member))),
                      time = factor(seq_len(n_time)),
                      group = factor(seq_len(2)))
    dd$y <- 0

    K <- matrix(c(1.00, 0.25, 0.10,
                  0.25, 1.44, 0.20,
                  0.10, 0.20, 0.81), n_member, n_member)
    dimnames(K) <- list(levels(dd$member), levels(dd$member))
    phi <- 0.45
    extra_sd <- 1.3
    R_member <- cov2cor(K)
    R_time <- outer(seq_len(n_time), seq_len(n_time),
                    function(i, j) phi^abs(i - j))
    sd_member <- unname(sqrt(diag(K)) * extra_sd)
    R_full <- kronecker(R_time, R_member)
    sd_full <- rep(sd_member, n_time)

    form <- y ~ 1 +
        separable(propto(0 + member, K) %x% ar1(0 + time) | group)
    dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
    env <- list2env(list(K = K), parent = environment())
    environment(form) <- env
    environment(dense_form) <- env
    case <- list(
        form = form,
        dense_form = dense_form,
        dd = dd,
        theta = c(log(extra_sd), ar1_to_theta(phi)),
        theta_dense = c(log(sd_full), put_cor(R_full)),
        R_full = R_full,
        sd_full = sd_full,
        codes = unname(c(.valid_covstruct[["propto"]],
                         .valid_covstruct[["ar1"]])),
        kinds = c(sep_kind_code("propto"), sep_kind_code("ar1")),
        scale_kinds = c(1L, 1L),
        dispatch = 1L,
        scale_mode = 1L,
        scale_spec = 0L
    )

    expect_separable_case_vc(case)
    expect_separable_case_nll(case)
})

test_that("separable supports equalto matrix margins without estimated scale", {
    n_member <- 3
    n_time <- 4
    dd <- expand.grid(member = factor(paste0("m", seq_len(n_member))),
                      time = factor(seq_len(n_time)),
                      group = factor(seq_len(2)))
    dd$y <- 0

    K <- matrix(c(1.00, 0.25, 0.10,
                  0.25, 1.44, 0.20,
                  0.10, 0.20, 0.81), n_member, n_member)
    phi <- 0.35
    R_member <- cov2cor(K)
    R_time <- outer(seq_len(n_time), seq_len(n_time),
                    function(i, j) phi^abs(i - j))
    sd_member <- unname(sqrt(diag(K)))
    R_full <- kronecker(R_time, R_member)
    sd_full <- rep(sd_member, n_time)

    form <- y ~ 1 +
        separable(equalto(0 + member, K) %x% ar1(0 + time) | group)
    dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
    env <- list2env(list(K = K), parent = environment())
    environment(form) <- env
    environment(dense_form) <- env
    case <- list(
        form = form,
        dense_form = dense_form,
        dd = dd,
        theta = ar1_to_theta(phi),
        theta_dense = c(log(sd_full), put_cor(R_full)),
        R_full = R_full,
        sd_full = sd_full,
        codes = unname(c(.valid_covstruct[["equalto"]],
                         .valid_covstruct[["ar1"]])),
        kinds = c(sep_kind_code("equalto"), sep_kind_code("ar1")),
        scale_kinds = c(0L, 1L),
        dispatch = 1L,
        scale_mode = 0L,
        scale_spec = integer()
    )

    expect_separable_case_vc(case)
    expect_separable_case_nll(case)
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

    case <- make_sep_margin_chain_case(
        list(make_dense_margin("us", n = 2),
             make_ar1_margin("hetar1", n = 3),
             make_toep_margin("homtoep", n = 2)),
        vars = c("member", "time", "item"),
        scale_mode = "product"
    )
    fit3 <- fit_fixed_theta(case$form, case$dd, case$theta)
    sims3 <- simulate(fit3, nsim = 1)
    expect_s3_class(sims3, "data.frame")
    expect_equal(dim(sims3), c(nrow(case$dd), 1L))
    expect_true(is.numeric(sims3[[1]]))
})
