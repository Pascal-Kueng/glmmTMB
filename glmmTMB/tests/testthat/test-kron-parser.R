stopifnot(require("testthat"), require("glmmTMB"))

test_that("kron parser preserves the product specification", {
    split <- reformulas::splitForm(
        y ~ kron((us(0 + member) %x% ar1(0 + time)) %x%
                     diag(0 + outcome) | site),
        specials = "kron"
    )
    expect_identical(split$reTrmClasses, "kron")
    expect_identical(
        split$reTrmFormulas[[1L]],
        quote((us(0 + member) %x% ar1(0 + time)) %x%
                  diag(0 + outcome) | site)
    )
    expect_identical(split$reTrmAddArgs[[1L]], quote(kron()))

    spec <- glmmTMB:::.kron_parse(
        split$reTrmFormulas[[1L]],
        split$reTrmAddArgs[[1L]]
    )

    expect_identical(spec$struc, c("us", "ar1", "diag"))
    expect_identical(
        spec$expr,
        list(quote(0 + member), quote(0 + time), quote(0 + outcome))
    )
})

test_that("kron parser rejects malformed product terms", {
    expect_error(
        glmmTMB:::.kron_parse(
            quote(us(0 + member) %x% ar1(0 + time)),
            quote(kron())
        ),
        "must contain a product random-effects term"
    )
    expect_error(
        glmmTMB:::.kron_parse(
            quote(us(0 + member) | group),
            quote(kron())
        ),
        "needs at least two margins"
    )
    expect_error(
        glmmTMB:::.kron_parse(
            quote(toep(0 + member) %x% ar1(0 + time) | group),
            quote(kron())
        ),
        "unsupported kron\\(\\) margin: toep"
    )
    expect_error(
        glmmTMB:::.kron_parse(
            quote(us(x = 0 + member) %x% ar1(0 + time) | group),
            quote(kron())
        ),
        "margins must look like"
    )
    expect_error(
        glmmTMB:::.kron_parse(
            quote(us(0 + member) %x% ar1(0 + time) | group),
            quote(kron(scale = "global"))
        ),
        "does not accept additional arguments"
    )
})
