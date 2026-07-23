## Workaround to associate numeric values with factor levels in a way
## that survives through the lme4 machinery.

##' Create a factor with numeric interpretable factor levels.
##'
##' Some \code{glmmTMB} covariance structures require extra
##' information, such as temporal or spatial
##' coordinates. \code{numFactor} allows to associate such extra
##' information as part of a factor via the factor levels. The
##' original numeric coordinates are recoverable without loss of
##' precision using the function \code{parseNumLevels}.  Factor levels
##' are sorted coordinate wise from left to right: first coordinate is
##' fastest running.
##' @title Factor with numeric interpretable levels.
##' @param x Vector, matrix or data.frame that constitute the
##'     coordinates.
##' @param ... Additional vectors, matrices or data.frames that
##'     constitute the coordinates.
##' @return Factor with specialized coding of levels.
##' @examples
##' ## 1D example
##' numFactor(sample(1:5,20,TRUE))
##' ## 2D example
##' coords <- cbind( sample(1:5,20,TRUE), sample(1:5,20,TRUE) )
##' (f <- numFactor(coords))
##' parseNumLevels(levels(f)) ## Sorted
##' ## Used as part of a model.matrix
##' model.matrix( ~f )
##' ## parseNumLevels( colnames(model.matrix( ~f )) )
##' ## Error: 'Failed to parse numeric levels: (Intercept)'
##' parseNumLevels( colnames(model.matrix( ~ f-1 )) )
##' @export
numFactor <- function(x, ...) {
    y <- data.frame(x, ...)
    if( !all( sapply(y, is.numeric) | sapply(y, is.factor)) )
        stop("All arguments to 'numFactor' must be numeric or factor.")
    asChar <- function(y) {
        y <- lapply(y, as.character)
        ans <- do.call("paste", c(y, list(sep=",")))
        paste0("(", ans, ")")
    }
    fac <- asChar(y)
    ndup <- !duplicated(fac)
    y0 <- y[ndup, , drop=FALSE]
    for (col in seq_along(y0) ) {
        y0 <- y0[ order( y0[[col]] ), , drop=FALSE]
    }
    facLevels <- asChar(y0)
    factor( fac, levels = facLevels )
}

##' @rdname numFactor
##' @param levels Character vector to parse into numeric values.
##' @importFrom stats complete.cases
##' @export
parseNumLevels <- function(levels) {
    ## Strip initial (irrelevant) characters:
    tmp <- sub("^.*(\\(.+\\))$", "\\1", levels)
    ## Now tmp must have the form ([0-9]*,[0-9]*,...)
    ## Otherwise it's an error
    tmp <- sub("^\\(", "", tmp)
    tmp <- sub("\\)$", "", tmp)
    ## Split string and convert to numeric
    ans <- lapply( strsplit(tmp, ","), as.numeric )
    ans <- t( do.call("cbind", ans) )
    ## if(any(is.na(ans))) stop("Failed to parse numeric levels.")
    if(any(is.na(ans))) {
        stop("Failed to parse numeric levels: ",
             levels[!complete.cases(ans)])
    }
    ans
}

.kron_supported <- c("homdiag", "diag", "homcs", "cs", "us", "ar1",
                     "hetar1")

.kron_call_name <- function(x) {
    if (!is.call(x)) return("")
    deparse1(x[[1L]], collapse = "")
}

.kron_flatten <- function(x) {
    while (is.call(x) && identical(.kron_call_name(x), "(") &&
           length(x) == 2L) x <- x[[2L]]
    if (is.call(x) && identical(.kron_call_name(x), "%x%") &&
        length(x) == 3L) {
        return(c(.kron_flatten(x[[2L]]), .kron_flatten(x[[3L]])))
    }
    list(x)
}

.kron_shape_npar <- function(struc, dim) {
    if (dim < 1L) stop("kron() margins cannot have zero columns",
                       call. = FALSE)
    if (dim < 2L && struc %in% c("homcs", "cs", "ar1", "hetar1")) {
        stop(struc, "() needs at least two margin columns",
             call. = FALSE)
    }
    as.integer(switch(struc,
                      homdiag = 0L,
                      diag = dim - 1L,
                      homcs = 1L,
                      cs = dim,
                      us = dim - 1L + dim * (dim - 1L) / 2L,
                      ar1 = 1L,
                      hetar1 = dim,
                      stop("unsupported kron() margin: ", struc,
                           call. = FALSE)))
}

.kron_parse <- function(bar, add_arg) {
    if (!is.call(bar) || !identical(.kron_call_name(bar), "|") ||
        length(bar) != 3L) {
        stop("kron() must contain a product random-effects term",
             call. = FALSE)
    }
    if (!is.call(add_arg) || !identical(.kron_call_name(add_arg), "kron") ||
        length(add_arg) != 1L) {
        stop("kron() does not accept additional arguments", call. = FALSE)
    }

    margins <- .kron_flatten(bar[[2L]])
    if (length(margins) < 2L) {
        stop("kron() needs at least two margins joined by %x%",
             call. = FALSE)
    }
    parsed <- lapply(margins, function(x) {
        named <- !is.null(names(x)) && any(nzchar(names(x)[-1L]))
        if (!is.call(x) || length(x) != 2L || named) {
            stop("kron() margins must look like us(1 + x) or ",
                 "ar1(0 + time)", call. = FALSE)
        }
        struc <- .kron_call_name(x)
        if (!struc %in% .kron_supported) {
            stop("unsupported kron() margin: ", struc, call. = FALSE)
        }
        list(struc = struc, expr = x[[2L]])
    })
    expr <- lapply(parsed, `[[`, "expr")
    margin_names <- vapply(expr, deparse1, character(1), collapse = "")
    list(struc = vapply(parsed, `[[`, character(1), "struc"),
         expr = expr, margin_names = margin_names)
}

.kron_frame_formula <- function(x) {
    if (is.name(x) || !is.language(x)) return(x)
    if (is.call(x) && identical(.kron_call_name(x), "kron")) {
        if (length(x) != 2L) {
            stop("kron() accepts exactly one product random-effects term",
                 call. = FALSE)
        }
        bar <- x[[2L]]
        spec <- .kron_parse(bar, quote(kron()))
        return(Reduce(function(a, b) call("+", a, b),
                      c(spec$expr, list(bar[[3L]]))))
    }
    for (i in seq_along(x)[-1L]) x[[i]] <- .kron_frame_formula(x[[i]])
    x
}

.kron_margin_matrix <- function(expr, fr, env) {
    form <- stats::as.formula(as.call(list(as.name("~"), expr)), env = env)
    X <- tryCatch(
        Matrix::sparse.model.matrix(form, data = fr),
        error = function(e) Matrix::Matrix(
            stats::model.matrix(form, data = fr), sparse = TRUE
        )
    )
    if (!ncol(X)) stop("kron() margins cannot have zero columns",
                       call. = FALSE)
    X
}

.kron_product_names <- function(x) {
    Reduce(function(a, b) as.vector(outer(a, b, paste, sep = ":")), x)
}

.kron_build <- function(spec, fr, env) {
    X <- lapply(spec$expr, .kron_margin_matrix, fr = fr, env = env)
    dims <- vapply(X, ncol, integer(1))
    structured <- spec$struc %in% c("ar1", "hetar1")
    one_hot <- vapply(X, function(x) {
        all(Matrix::rowSums(x != 0) <= 1L) && all(x@x == 1)
    }, logical(1))
    if (any(structured & !one_hot)) {
        bad <- which(structured & !one_hot)[[1L]]
        stop(spec$struc[[bad]], "() needs a no-intercept indicator margin",
             call. = FALSE)
    }

    product_t <- Reduce(
        function(ans, x) Matrix::KhatriRao(Matrix::t(x), ans),
        X[-1L], init = Matrix::t(X[[1L]])
    )
    margin_columns <- lapply(X, colnames)
    cnms <- .kron_product_names(margin_columns)
    rownames(product_t) <- cnms
    npar <- vapply(seq_along(dims), function(i) {
        .kron_shape_npar(spec$struc[[i]], dims[[i]])
    }, integer(1))
    source_vars <- unique(unlist(lapply(spec$expr, all.vars),
                                 use.names = FALSE))
    source_vars <- source_vars[vapply(source_vars, function(nm) {
        value <- if (nm %in% names(fr)) fr[[nm]] else tryCatch(
            eval(as.name(nm), envir = env), error = function(e) NULL
        )
        !is.null(value) && NROW(value) == nrow(fr)
    }, logical(1))]
    info <- structure(
        list(kronDims = as.integer(dims),
             kronCodes = as.integer(unname(.valid_covstruct[spec$struc])),
             kronMarginNames = spec$margin_names,
             kronMarginColumns = margin_columns,
             kronSourceVars = source_vars,
             ntheta = as.integer(1L + sum(npar))),
        class = "glmmTMB_kron_spec"
    )
    list(product_t = product_t, cnms = cnms, info = info)
}

.kron_calls <- function(x) {
    if (is.name(x) || !is.language(x)) return(list())
    ans <- if (is.call(x) && identical(.kron_call_name(x), "kron")) {
        list(x)
    } else {
        list()
    }
    children <- lapply(as.list(x)[-1L], .kron_calls)
    c(ans, unlist(children, recursive = FALSE))
}

.kron_formula_info <- function(forms) {
    calls <- unlist(lapply(forms, .kron_calls), recursive = FALSE)
    if (!length(calls)) {
        return(list(has_kron = FALSE, source_vars = character(),
                    margin_vars = character(),
                    frame_vars = character()))
    }
    specs <- lapply(calls, function(x) {
        if (length(x) != 2L) {
            stop("kron() accepts exactly one product random-effects term",
                 call. = FALSE)
        }
        .kron_parse(x[[2L]], quote(kron()))
    })
    expr <- unlist(lapply(specs, `[[`, "expr"), recursive = FALSE)
    margin_vars <- unique(unlist(lapply(expr, all.vars), use.names = FALSE))
    frame_vars <- unlist(lapply(expr, function(x) {
        form <- stats::as.formula(as.call(list(as.name("~"), x)))
        vars <- as.list(attr(stats::terms(form), "variables"))[-1L]
        vapply(vars, function(v) {
            if (is.name(v)) as.character(v) else deparse1(v, collapse = "")
        }, character(1))
    }), use.names = FALSE)
    list(has_kron = TRUE,
         source_vars = unique(unlist(lapply(
             calls, function(x) all.vars(x[[2L]])
         ), use.names = FALSE)),
         margin_vars = margin_vars,
         frame_vars = unique(frame_vars))
}
