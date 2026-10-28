## Tests for the RTMB kron() (Kronecker product) covariance structure

context("RTMB kron() covariance structure")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-6

set.seed(6001)
## two observations per member, time and dyad, so that the random effects are
## not confounded with the residual
kron_data <- expand.grid(
  member = factor(1:2),
  time = factor(1:5),
  dyad = factor(1:25),
  obs = 1:2
)
S_member <- matrix(c(1, 0.5, 0.5, 1.5), 2)
S_time <- 0.6^abs(outer(1:5, 1:5, "-"))
## member varies fastest in kron_data's rows, so S_member comes last
L <- t(chol(kronecker(S_time, S_member)))
b <- as.vector(L %*% matrix(rnorm(10 * 25), 10))
kron_data$x <- rnorm(nrow(kron_data))
kron_data$y <- 1 + 0.5 * kron_data$x + rep(b, 2) + rnorm(nrow(kron_data), sd = 0.5)

fit_rtmb <- function(formula, data = kron_data) {
  glmmTMB(formula, data = data, control = glmmTMBControl(use_rtmb = TRUE))
}

m_kron <- fit_rtmb(y ~ x + kron(us(0 + member) %x% ar1(0 + time) | dyad))

test_that("kron() matches an equivalent ordinary structure", {
  ## homdiag() as a later margin is the identity (its SD is fixed at 1)
  f1 <- fit_rtmb(y ~ kron(us(0 + member) %x% homdiag(0 + time) | dyad))
  f2 <- fit_rtmb(y ~ us(0 + member | dyad:time))
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
  expect_equal(unname(getME(f1, "theta")), unname(getME(f2, "theta")),
               tolerance = 1e-4)
})

test_that("kron() Gaussian likelihood matches the dense marginal likelihood", {
  vc <- VarCorr(m_kron)$cond
  expect_equal(names(vc), c("dyad (kron 1/2)", "dyad (kron 2/2)"))
  ## covariance of one dyad's random effects, in Z's column order: time (the
  ## last margin) varies fastest
  S_block <- kronecker(vc[[1]], vc[[2]])
  Z <- as.matrix(getME(m_kron, "Z"))
  V <- Z %*% kronecker(diag(25), S_block) %*% t(Z) +
    sigma(m_kron)^2 * diag(nrow(kron_data))
  res <- kron_data$y - getME(m_kron, "X") %*% fixef(m_kron)$cond
  R <- chol(V)
  z <- backsolve(R, res, transpose = TRUE) ## whitened residuals, iid N(0, 1)
  nll <- sum(log(diag(R))) - sum(dnorm(z, log = TRUE))
  expect_equal(-as.numeric(logLik(m_kron)), nll, tolerance = tol_logLik)
})

test_that("three-margin kron() matches ar1()", {
  dd <- expand.grid(a = factor(1:2), time = factor(1:4), b = factor(1:3),
                    g = factor(1:10))
  ## rows of dd: a varies fastest, then time, then b
  S <- kronecker(diag(3), kronecker(0.5^abs(outer(1:4, 1:4, "-")), diag(2)))
  dd$y <- as.vector(t(chol(S)) %*% matrix(rnorm(24 * 10), 24)) +
    rnorm(nrow(dd), sd = 0.5)
  f1 <- fit_rtmb(y ~ kron(homdiag(0 + a) %x% ar1(0 + time) %x% homdiag(0 + b) | g), dd)
  f2 <- fit_rtmb(y ~ ar1(0 + time | g:a:b), dd)
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
})

test_that("kron() margins can have intercepts and slopes", {
  ## homdiag() as a later margin is the identity (its SD is fixed at 1)
  f1 <- fit_rtmb(y ~ x + kron(us(1 + x) %x% homdiag(0 + time) | dyad))
  f2 <- fit_rtmb(y ~ x + us(1 + x | dyad:time))
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
})

test_that("kron() margins can be propto or equalto", {
  ## an equalto first margin has no free scale, so the next margin keeps its own
  f1 <- fit_rtmb(y ~ kron(equalto(0 + member, S_member) %x% ar1(0 + time) | dyad))
  expect_equal(c(VarCorr(f1)$cond[[1]]), c(S_member))
  expect_length(getME(f1, "theta"), 2)
  ## a kron() term may have no estimated parameters, also when a term follows it
  f2 <- fit_rtmb(y ~ kron(equalto(0 + member, S_member) %x% equalto(0 + time, S_time) | dyad) +
                   (1 | dyad))
  expect_length(getME(f2, "theta"), 1)
  ## from glmmTMB#1353: diag(0 + DOP) is a 1 x 1 slope margin, so this is
  ## propto(0 + spp:DOP | site:sample, C)
  data("Salamanders", package = "glmmTMB")
  lev <- levels(Salamanders$spp)
  C_spp <- 0.5^abs(outer(seq_along(lev), seq_along(lev), "-"))
  dimnames(C_spp) <- list(lev, lev)
  C_sppDOP <- C_spp
  dimnames(C_sppDOP) <- rep(list(paste0("spp", lev, ":DOP")), 2)
  ## glmmTMB() looks up the matrices from where it is called, so call it here
  ctl <- glmmTMBControl(use_rtmb = TRUE)
  f3 <- glmmTMB(count ~ spp + kron(diag(0 + DOP) %x% propto(0 + spp, C_spp) | site:sample),
                data = Salamanders, family = poisson, control = ctl)
  f4 <- glmmTMB(count ~ spp + propto(0 + spp:DOP | site:sample, C_sppDOP),
                data = Salamanders, family = poisson, control = ctl)
  expect_equal(logLik(f3), logLik(f4), tolerance = tol_logLik)
  ## as a later margin, propto's scale is fixed
  expect_length(getME(f3, "theta"), 1)
})

test_that("kron() margins can be rr", {
  f1 <- fit_rtmb(y ~ kron(rr(0 + member, 1) %x% homdiag(0 + time) | dyad))
  f2 <- fit_rtmb(y ~ rr(0 + member | dyad:time, 1))
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
})

test_that("kron() works next to other terms", {
  dd <- transform(kron_data, site = factor(as.integer(dyad) %% 5))
  dd$y <- dd$y + rnorm(5)[dd$site]
  f1 <- fit_rtmb(y ~ x + (1 | site) + kron(homdiag(0 + member) %x% ar1(0 + time) | dyad), dd)
  f2 <- fit_rtmb(y ~ x + (1 | site) + ar1(0 + time | dyad:member), dd)
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
  expect_equal(names(VarCorr(f1)$cond), c("site", "dyad (kron 1/2)", "dyad (kron 2/2)"))
})

test_that("kron() expands nested grouping factors", {
  dd <- transform(kron_data, site = factor(as.integer(dyad) %% 5))
  nested <- glmmTMB(
    y ~ x + kron(us(0 + member) %x% ar1(0 + time) | site/dyad),
    data = dd, doFit = FALSE)
  explicit <- glmmTMB(
    y ~ x + kron(us(0 + member) %x% ar1(0 + time) | dyad:site) +
      kron(us(0 + member) %x% ar1(0 + time) | site),
    data = dd, doFit = FALSE)
  expect_equal(nested$condReStruc, explicit$condReStruc)
  expect_equal(nested$data.tmb$Z, explicit$data.tmb$Z)
})

test_that("kron() works with formula(), simulate() and predict(newdata)", {
  expect_identical(deparse1(formula(m_kron)),
                   "y ~ x + kron(us(0 + member) %x% ar1(0 + time) | dyad)")
  s <- simulate(m_kron, nsim = 2, seed = 1)
  expect_equal(dim(s), c(nrow(kron_data), 2L))
  expect_true(all(is.finite(unlist(s))))
  expect_equal(predict(m_kron, newdata = kron_data[1:10, ]),
               predict(m_kron)[1:10])

  obj <- m_kron$obj
  set_simcodes(obj, "zero")
  expect_equal(obj$simulate()$b, rep(0, 10 * 25))
  set_simcodes(obj, "fix")
  expect_equal(obj$simulate()$b, obj$report()$b)
  set_simcodes(obj, "random") ## restore the default
})

test_that("kron() gives clear errors", {
  expect_error(
    glmmTMB(y ~ kron(us(0 + member) %x% ar1(0 + time) | dyad),
            data = kron_data, control = glmmTMBControl(use_rtmb = FALSE)),
    "needs the RTMB back-end")
  expect_error(fit_rtmb(y ~ kron(us(0 + member) %x% kron(0 + time) | dyad)),
               "cannot be kron")
  expect_error(fit_rtmb(y ~ kron(us(0 + member) %x% time | dyad)),
               "must look like")
  expect_error(fit_rtmb(y ~ kron(us(0 + member | dyad) %x% ar1(0 + time | dyad) | dyad)),
               "no grouping bar")
  expect_error(fit_rtmb(y ~ kron(us(0 + member) %x% ar1(0 + time))),
               "needs a grouping factor")
  expect_error(fit_rtmb(y ~ kron(us(0 + member) %x% ar1(0 + time) | 1)),
               "factor with one level")
})

test_that("kron() rejects unsupported Kenward-Roger inference", {
  m_reml <- update(m_kron, data = kron_data, REML = TRUE)
  expect_error(dof_KR(m_reml), "not supported.*kron")
  expect_error(summary(m_reml, ddf = "kenward-roger"), "not supported.*kron")
})
