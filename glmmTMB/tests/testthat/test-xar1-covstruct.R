stopifnot(require("testthat"),
          require("glmmTMB"))

mk_xar1_dat <- function(ng = 3, nt = 3, nmember = 2, drop = FALSE) {
    dd <- expand.grid(group = factor(seq_len(ng)),
                      member = seq_len(nmember),
                      time = seq_len(nt))
    if (drop) {
        dd <- dd[!(dd$group == 1 & dd$member == 2 & dd$time == nt), ]
        full <- numFactor(rep(seq_len(nmember), nt),
                          rep(seq_len(nt), each = nmember))
        obs <- numFactor(dd$member, dd$time)
        dd$mt <- factor(as.character(obs), levels = levels(full))
    } else {
        dd$mt <- numFactor(dd$member, dd$time)
    }
    set.seed(101)
    dd$y <- rnorm(nrow(dd))
    dd
}

rho_to_theta <- function(rho) rho / sqrt(1 - rho^2)
homcs_to_theta <- function(rho, n) {
    a <- 1 / (n - 1)
    qlogis((rho + a) / (1 + a))
}

xar1_report <- function(obj) {
    tmb <- fitTMB(obj, doOptim = FALSE)
    tmb$env$report(tmb$env$last.par.best)
}

test_that("separable AR1 structures parse and evaluate", {
    dd <- mk_xar1_dat()

    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0, doFit = FALSE))
    expect_equal(m1$condReStruc[[1]]$blockNumTheta, 3)
    expect_equal(m1$condReStruc[[1]]$sepMembers, c(1, 2, 1, 2, 1, 2))
    expect_equal(m1$condReStruc[[1]]$sepTimes, c(1, 1, 2, 2, 3, 3))
    obj1 <- fitTMB(m1, doOptim = FALSE)
    expect_true(is.finite(obj1$fn(obj1$par)))

    m2 <- glmmTMB(y ~ 1 + unxar1(mt + 0 | group), data = dd,
                  dispformula = ~0, doFit = FALSE)
    expect_equal(m2$condReStruc[[1]]$blockNumTheta, 4)
    obj2 <- fitTMB(m2, doOptim = FALSE)
    expect_true(is.finite(obj2$fn(obj2$par)))
})

test_that("separable AR1 structures support more than two member levels", {
    dd <- mk_xar1_dat(nmember = 3)

    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0, doFit = FALSE))
    expect_equal(m1$condReStruc[[1]]$blockSize, 9)
    expect_equal(m1$condReStruc[[1]]$blockNumTheta, 3)
    expect_equal(m1$condReStruc[[1]]$sepMembers, rep(1:3, 3))
    expect_equal(m1$condReStruc[[1]]$sepTimes, rep(1:3, each = 3))
    obj1 <- fitTMB(m1, doOptim = FALSE)
    expect_true(is.finite(obj1$fn(obj1$par)))

    m2 <- glmmTMB(y ~ 1 + unxar1(mt + 0 | group), data = dd,
                  dispformula = ~0, doFit = FALSE)
    expect_equal(m2$condReStruc[[1]]$blockSize, 9)
    expect_equal(m2$condReStruc[[1]]$blockNumTheta, 7)
    obj2 <- fitTMB(m2, doOptim = FALSE)
    expect_true(is.finite(obj2$fn(obj2$par)))
})

test_that("membertime and mt helpers provide inline member-time coordinates", {
    dd <- mk_xar1_dat()
    dd$role <- factor(dd$member, labels = c("A", "B"))
    dd$day <- dd$time

    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(membertime(role, day) + 0 | group),
                data = dd, dispformula = ~0, doFit = FALSE))
    expect_equal(m1$condReStruc[[1]]$blockSize, 6)
    expect_equal(m1$condReStruc[[1]]$sepMembers, c(1, 2, 1, 2, 1, 2))
    expect_equal(m1$condReStruc[[1]]$sepTimes, c(1, 1, 2, 2, 3, 3))

    m2 <- glmmTMB(y ~ 1 + unxar1(mt(role, day) + 0 | group),
                  data = dd, dispformula = ~0, doFit = FALSE)
    expect_equal(m2$condReStruc[[1]]$blockNumTheta, 4)
})

test_that("separable AR1 covariance matrices have expected pattern", {
    dd <- mk_xar1_dat()
    phi <- 0.5
    rho <- 0.25

    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0, doFit = FALSE,
                start = list(theta = c(log(2), homcs_to_theta(rho, 2),
                                       rho_to_theta(phi)))))
    r1 <- xar1_report(m1)
    cc1 <- r1$corr[[1]]
    expect_equal(cc1[1, 3], phi, tolerance = 1e-8)
    expect_equal(cc1[1, 2], rho, tolerance = 1e-8)
    expect_equal(cc1[1, 4], rho * phi, tolerance = 1e-8)
    expect_equal(r1$sd[[1]], rep(2, 6), tolerance = 1e-8)

    m2 <- glmmTMB(y ~ 1 + unxar1(mt + 0 | group), data = dd,
                  dispformula = ~0, doFit = FALSE,
                  start = list(theta = c(log(2), log(3), rho_to_theta(rho),
                                         rho_to_theta(phi))))
    r2 <- xar1_report(m2)
    cc2 <- r2$corr[[1]]
    expect_equal(cc2[1, 3], phi, tolerance = 1e-8)
    expect_equal(cc2[1, 2], rho, tolerance = 1e-8)
    expect_equal(cc2[1, 4], rho * phi, tolerance = 1e-8)
    expect_equal(r2$sd[[1]], rep(c(2, 3), 3), tolerance = 1e-8)
})

test_that("separable AR1 covariance matrices reduce correctly at zero correlations", {
    dd <- mk_xar1_dat()
    make_hom <- function(rho, phi) {
        suppressWarnings(
            glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                    dispformula = ~0, doFit = FALSE,
                    start = list(theta = c(log(1), homcs_to_theta(rho, 2),
                                           rho_to_theta(phi)))))
    }

    cc_rho0 <- xar1_report(make_hom(rho = 0, phi = 0.4))$corr[[1]]
    expect_equal(cc_rho0[1, 2], 0, tolerance = 1e-8)
    expect_equal(cc_rho0[1, 4], 0, tolerance = 1e-8)
    expect_equal(cc_rho0[1, 3], 0.4, tolerance = 1e-8)

    cc_phi0 <- xar1_report(make_hom(rho = 0.3, phi = 0))$corr[[1]]
    expect_equal(cc_phi0[1, 2], 0.3, tolerance = 1e-8)
    expect_equal(cc_phi0[1, 3], 0, tolerance = 1e-8)
    expect_equal(cc_phi0[1, 4], 0, tolerance = 1e-8)

    cc_both0 <- xar1_report(make_hom(rho = 0, phi = 0))$corr[[1]]
    expect_equal(cc_both0[lower.tri(cc_both0)], rep(0, 15), tolerance = 1e-8)
})

test_that("unxar1 with equal member SDs matches homcsxar1", {
    dd <- mk_xar1_dat()
    rho <- 0.35
    phi <- -0.25
    sd <- 1.4

    m_hom <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0, doFit = FALSE,
                start = list(theta = c(log(sd), homcs_to_theta(rho, 2),
                                       rho_to_theta(phi)))))
    m_un <- glmmTMB(y ~ 1 + unxar1(mt + 0 | group), data = dd,
                    dispformula = ~0, doFit = FALSE,
                    start = list(theta = c(log(sd), log(sd), rho_to_theta(rho),
                                           rho_to_theta(phi))))

    r_hom <- xar1_report(m_hom)
    r_un <- xar1_report(m_un)
    expect_equal(r_hom$sd[[1]], r_un$sd[[1]], tolerance = 1e-8)
    expect_equal(r_hom$corr[[1]], r_un$corr[[1]], tolerance = 1e-8)
})

test_that("separable AR1 validates coordinate input", {
    dd <- mk_xar1_dat()
    expect_error(
        glmmTMB(y ~ 1 + homcsxar1(mt | group), data = dd, doFit = FALSE),
        "without an intercept"
    )

    dd$bad1 <- numFactor(dd$time)
    expect_error(
        glmmTMB(y ~ 1 + homcsxar1(bad1 + 0 | group), data = dd, doFit = FALSE),
        "exactly two coordinates"
    )

    dd$badtime <- numFactor(dd$member, dd$time * 2)
    expect_error(
        glmmTMB(y ~ 1 + homcsxar1(badtime + 0 | group), data = dd, doFit = FALSE),
        "discrete unit-spaced time"
    )

    dd_dup <- dd
    dd_dup$role <- factor(ifelse(dd_dup$member == 1, "male", "female"),
                          levels = c("male", "female"))
    dd_dup$role[dd_dup$group == 1 & dd_dup$time == 1] <- "male"
    expect_error(
        glmmTMB(y ~ 1 + unxar1(membertime(role, time) + 0 | group),
                data = dd_dup, dispformula = ~0, doFit = FALSE),
        "duplicate observations.*same-sex dyads cannot use gender"
    )
})

test_that("separable AR1 handles missing observations within groups", {
    dd <- mk_xar1_dat(drop = TRUE)
    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0, doFit = FALSE))
    expect_equal(m1$condReStruc[[1]]$blockSize, 6)
    obj1 <- fitTMB(m1, doOptim = FALSE)
    expect_true(is.finite(obj1$fn(obj1$par)))

    dd_unit <- mk_xar1_dat()
    full <- numFactor(rep(1:2, 3), rep(1:3, each = 2))
    dd_unit <- dd_unit[!(dd_unit$group == 1 & dd_unit$member == 2), ]
    obs <- numFactor(dd_unit$member, dd_unit$time)
    dd_unit$mt <- factor(as.character(obs), levels = levels(full))
    expect_warning(
        glmmTMB(y ~ 1 + unxar1(mt + 0 | group), data = dd_unit,
                dispformula = ~0, doFit = FALSE),
        "fewer than two member/role levels"
    )
})

test_that("separable AR1 warns about dispersion interpretation", {
    dd <- mk_xar1_dat()
    expect_warning(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                doFit = FALSE),
        "independent residual/nugget variance"
    )
    expect_warning(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~ member, doFit = FALSE),
        "structured dispersion model"
    )
    expect_warning(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0, doFit = FALSE),
        "same individual within each group across time"
    )

    dd$count <- rpois(nrow(dd), 2)
    expect_warning(
        glmmTMB(count ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                family = poisson, dispformula = ~ member, doFit = FALSE),
        "dispformula is ignored"
    )

    dd$count_nb <- rnbinom(nrow(dd), mu = 2, size = 5)
    expect_warning(
        glmmTMB(count_nb ~ 1 + unxar1(mt + 0 | group), data = dd,
                family = nbinom2, dispformula = ~ member, doFit = FALSE),
        "family-specific dispersion parameter"
    )
})

test_that("xar1 structures print compactly for model and summary", {
    dd <- mk_xar1_dat(ng = 20)
    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0))
    p1 <- capture.output(print(m1))
    s1 <- capture.output(print(summary(m1)))
    expect_true(any(grepl("mt member", p1, fixed = TRUE)))
    expect_true(any(grepl("time lag 1", p1, fixed = TRUE)))
    expect_true(any(grepl("(xar1)", p1, fixed = TRUE)))
    expect_true(any(grepl("mt member", s1, fixed = TRUE)))
    expect_true(any(grepl("time lag 1", s1, fixed = TRUE)))
    expect_false(any(grepl("mt\\(1,3\\)", p1)))
    expect_false(any(grepl("mt\\(1,3\\)", s1)))

    v1 <- capture.output(print(VarCorr(m1)))
    expect_true(any(grepl("mt\\(1,3\\)", v1)))
})

test_that("xar1 compact output combines with ordinary random effects", {
    dd <- mk_xar1_dat(ng = 20)
    dd$Idiff <- ifelse(dd$member == 1, -1, 1)
    dd$member_pos <- factor(dd$Idiff, levels = c(-1, 1),
                            labels = c("member_minus", "member_plus"))

    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + (1 | group) + (0 + Idiff | group) +
                    homcsxar1(mt(member_pos, time) + 0 | group),
                data = dd, dispformula = ~0))
    p1 <- capture.output(print(m1))
    s1 <- capture.output(print(summary(m1)))
    expect_true(any(grepl("(Intercept)", p1, fixed = TRUE)))
    expect_true(any(grepl("Idiff", p1, fixed = TRUE)))
    expect_true(any(grepl("member_pos x time member", p1, fixed = TRUE)))
    expect_true(any(grepl("member_pos x time member", s1, fixed = TRUE)))
})

test_that("inline membertime labels print as member x time", {
    dd <- mk_xar1_dat(ng = 20)
    dd$role <- factor(dd$member, labels = c("A", "B"))
    dd$day <- dd$time
    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(membertime(role, day) + 0 | group),
                data = dd, dispformula = ~0))
    p1 <- capture.output(print(m1))
    expect_true(any(grepl("role x day member", p1, fixed = TRUE)))
    expect_true(any(grepl("role x day time lag 1", p1, fixed = TRUE)))

    m2 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(mt(role, day) + 0 | group),
                data = dd, dispformula = ~0))
    p2 <- capture.output(print(m2))
    expect_true(any(grepl("role x day member", p2, fixed = TRUE)))

    m3 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(membertime(member = role, time = day) + 0 | group),
                data = dd, dispformula = ~0))
    p3 <- capture.output(print(m3))
    expect_true(any(grepl("role x day member", p3, fixed = TRUE)))
    expect_false(any(grepl("member = role x time = day", p3, fixed = TRUE)))
})

test_that("membertime factor levels are used as member labels in compact output", {
    dd <- mk_xar1_dat(ng = 20)
    dd$gender <- factor(dd$member, levels = 1:2, labels = c("female", "male"))
    dd$day <- dd$time

    m1 <- suppressWarnings(
        glmmTMB(y ~ 1 + homcsxar1(membertime(gender, day) + 0 | group),
                data = dd, dispformula = ~0))
    p1 <- capture.output(print(m1))
    expect_true(any(grepl("gender x day members female:male lag 1", p1, fixed = TRUE)))

    m2 <- glmmTMB(y ~ 1 + unxar1(membertime(gender, day) + 0 | group),
                  data = dd, dispformula = ~0)
    p2 <- capture.output(print(m2))
    s2 <- capture.output(print(summary(m2)))
    expect_true(any(grepl("gender x day female", p2, fixed = TRUE)))
    expect_true(any(grepl("gender x day male", p2, fixed = TRUE)))
    expect_true(any(grepl("gender x day members female:male lag 1", p2, fixed = TRUE)))
    expect_true(any(grepl("gender x day female", s2, fixed = TRUE)))
    expect_true(any(grepl("gender x day members female:male lag 1", s2, fixed = TRUE)))

    dd$gender_day <- membertime(dd$gender, dd$day)
    m3 <- glmmTMB(y ~ 1 + unxar1(gender_day + 0 | group),
                  data = dd, dispformula = ~0)
    p3 <- capture.output(print(m3))
    expect_true(any(grepl("gender_day female", p3, fixed = TRUE)))
    expect_true(any(grepl("gender_day members female:male lag 1", p3, fixed = TRUE)))
})

test_that("membertime numeric member labels preserve observed values", {
    dd <- mk_xar1_dat(ng = 20)
    dd$member01 <- dd$member - 1
    dd$member1020 <- ifelse(dd$member == 1, 10, 20)
    dd$day <- dd$time

    m1 <- glmmTMB(y ~ 1 + unxar1(membertime(member01, day) + 0 | group),
                  data = dd, dispformula = ~0)
    p1 <- capture.output(print(m1))
    expect_true(any(grepl("member01 x day 0", p1, fixed = TRUE)))
    expect_true(any(grepl("member01 x day 1", p1, fixed = TRUE)))
    expect_true(any(grepl("member01 x day members 0:1 lag 1", p1, fixed = TRUE)))

    m2 <- glmmTMB(y ~ 1 + unxar1(membertime(member1020, day) + 0 | group),
                  data = dd, dispformula = ~0)
    p2 <- capture.output(print(m2))
    expect_true(any(grepl("member1020 x day 10", p2, fixed = TRUE)))
    expect_true(any(grepl("member1020 x day 20", p2, fixed = TRUE)))
    expect_true(any(grepl("member1020 x day members 10:20 lag 1", p2, fixed = TRUE)))
})

test_that("dharma_xar1 exposes Gaussian marginal normalized residuals", {
    dd <- mk_xar1_dat(ng = 8, nt = 3)
    dd$x <- rnorm(nrow(dd))
    m1 <- suppressWarnings(
        glmmTMB(y ~ x + homcsxar1(mt + 0 | group), data = dd,
                dispformula = ~0))
    dx <- dharma_xar1(m1)

    expect_s3_class(dx, "glmmTMB")
    expect_false(is.null(attr(dx, "dharma_xar1")))
    expect_equal(attr(dx, "dharma_xar1")$type, "marginal")
    expect_equal(nobs(dx), nrow(dd))
    expect_equal(length(residuals(dx)), nobs(dx))
    expect_equal(length(predict(dx)), nobs(dx))
    expect_equal(predict(dx), as.vector(getME(m1, "X") %*% fixef(m1)$cond),
                 ignore_attr = TRUE)
    expect_gt(length(unique(predict(dx))), 1)
    expect_equal(nrow(model.frame(dx)), nobs(dx))

    sim <- simulate(dx, nsim = 3, seed = 1)
    expect_s3_class(sim, "data.frame")
    expect_equal(dim(sim), c(nobs(dx), 3))

    skip_if_not_installed("DHARMa")
    res <- DHARMa::simulateResiduals(dx, n = 5, seed = 1)
    expect_equal(length(res$scaledResiduals), nobs(dx))
    expect_equal(dim(res$simulatedResponse), c(nobs(dx), 5))
})

test_that("dharma_xar1 works for unxar1 and rejects non-Gaussian models", {
    dd <- mk_xar1_dat(ng = 6, nt = 3)
    m1 <- glmmTMB(y ~ 1 + unxar1(mt + 0 | group), data = dd,
                  dispformula = ~0)
    dx <- dharma_xar1(m1)
    expect_equal(nobs(dx), nrow(dd))
    expect_equal(length(residuals(dx)), nrow(dd))

    dd$count <- rpois(nrow(dd), lambda = 2)
    m2 <- suppressWarnings(
        glmmTMB(count ~ 1 + unxar1(mt + 0 | group), data = dd,
                family = poisson))
    expect_error(dharma_xar1(m2), "Gaussian identity-link")
})

test_that("membertime ignores unused factor levels when storing labels", {
    dd <- mk_xar1_dat(ng = 20)
    dd$role <- factor(dd$member, levels = c(1, 2, 3),
                      labels = c("female", "male", "unused"))
    dd$day <- dd$time

    mtval <- membertime(dd$role, dd$day)
    expect_equal(attr(mtval, "memberLabels"), c("1" = "female", "2" = "male"))

    m1 <- glmmTMB(y ~ 1 + unxar1(membertime(role, day) + 0 | group),
                  data = dd, dispformula = ~0)
    p1 <- capture.output(print(m1))
    expect_true(any(grepl("role x day female", p1, fixed = TRUE)))
    expect_false(any(grepl("unused", p1, fixed = TRUE)))
})
