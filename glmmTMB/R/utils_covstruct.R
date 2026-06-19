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
##'
##' \code{sepgrid} is similar to \code{numFactor}, but creates levels for
##' the complete Cartesian product of the supplied coordinate levels.  Factor
##' inputs preserve unused levels; non-factor inputs use sorted observed
##' values.  Use factors with explicit levels when globally unobserved cells
##' are part of the intended separable grid, for example an unobserved day in
##' an AR(1) time series.
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
##' @export
sepgrid <- function(x, ...) {
    ## Like numFactor(), but with levels for the complete Cartesian product.
    ## Factor inputs preserve unused levels; non-factors use sorted observed
    ## values.
    y <- data.frame(x, ...)

    ok <- vapply(y, function(z) is.numeric(z) || is.factor(z) ||
                   is.character(z) || is.integer(z), logical(1))
    if (!all(ok))
        stop("All arguments to 'sepgrid' must be numeric, factor, integer, or character.")

    ## Store coordinate indices so parseNumLevels() can recover dimensions.
    levs <- lapply(y, function(z) {
        if (is.factor(z)) levels(z) else sort(unique(z[!is.na(z)]))
    })
    vals <- Map(function(z, lev) match(if (is.factor(z)) as.character(z) else z, lev),
                y, levs)
    vals <- as.data.frame(vals)

    asChar <- function(y) {
        is_na <- !stats::complete.cases(y)
        y <- lapply(y, as.character)
        ans <- do.call("paste", c(y, list(sep=",")))
        ans <- paste0("(", ans, ")")
        ans[is_na] <- NA_character_
        ans
    }

    ## expand.grid() varies its first argument fastest; C++ uses the same order.
    grid <- do.call(expand.grid, c(lapply(levs, seq_along),
                                   list(KEEP.OUT.ATTRS = FALSE)))
    factor(asChar(vals), levels = asChar(grid))
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

## The helpers below compile the public product syntax
##
##   separable(us(0 + member) %x% ar1(0 + time) | group,
##             scale = us(0 + member))
##
## to an internal spec plus a lower-level random-effect term
##
##   separable(<combined margin variables> + 0 | group, <spec id>)
##
## `reformulas::splitForm()` still sees an ordinary covariance-structure
## special, but only carries a small spec id.  The structured margin/scale
## object is stored on the rewritten formula and passed explicitly to
## `getReStruc()`.  `mkReTrms()` initially sees the combined lower-level term,
## then glmmTMB replaces the separable term with the product of the marginal
## model matrices before the TMB structures are built.
.sep_deparse <- function(x) deparse1(x, collapse = "", width.cutoff = 500L)

.sep_call_name <- function(x) {
    if (!is.call(x)) return(NULL)
    .sep_deparse(x[[1]])
}

.sep_find_calls <- function(x, name) {
    ## Return all calls with the requested head.  This is used by glmmTMB() to
    ## identify generated `sepgrid(...)` model-frame columns that must keep their
    ## full Cartesian levels even when ordinary factors are level-dropped.
    if (!is.call(x)) return(list())
    ans <- if (identical(.sep_call_name(x), name)) list(x) else list()
    for (i in seq_along(x)[-1]) {
        ans <- c(ans, .sep_find_calls(x[[i]], name))
    }
    ans
}

.sepgrid_colnames <- function(...) {
    forms <- list(...)
    calls <- unlist(lapply(forms, function(f) {
        if (!inherits(f, "formula")) return(list())
        .sep_find_calls(f[[length(f)]], "sepgrid")
    }), recursive = FALSE)
    unique(vapply(calls, .sep_deparse, character(1)))
}

.sep_margin_varnames <- function(...) {
    forms <- list(...)
    vars <- unlist(lapply(forms, function(f) {
        specs <- attr(f, "separable_specs", exact = TRUE)
        if (is.null(specs)) return(character())
        unlist(lapply(specs, function(spec) {
            unique(unlist(lapply(spec$margins$expr, all.vars), use.names = FALSE))
        }), use.names = FALSE)
    }), use.names = FALSE)
    unique(vars)
}

.sep_is_zero <- function(x) {
    is.numeric(x) && length(x) == 1L && isTRUE(unname(x) == 0)
}

.sep_product_margin_label <- function(x) {
    if (is.call(x) && identical(.sep_call_name(x), "+") && length(x) == 3L) {
        if (.sep_is_zero(x[[2]]) && is.name(x[[3]])) return(.sep_deparse(x[[3]]))
        if (.sep_is_zero(x[[3]]) && is.name(x[[2]])) return(.sep_deparse(x[[2]]))
    }
    .sep_deparse(x)
}

.sep_product_margin_spec <- function(x) {
    if (!is.call(x) || length(x) != 2L)
        stop("separable() product margins must look like us(0 + role) ",
             "or ar1(0 + day).")
    data.frame(struc = .sep_deparse(x[[1]]),
               var = .sep_product_margin_label(x[[2]]),
               expr = I(list(x[[2]])),
               stringsAsFactors = FALSE)
}

.sep_add_calls <- function(x) {
    if (length(x) == 1L) return(x[[1]])
    Reduce(function(a, b) as.call(list(as.name("+"), a, b)), x)
}

.sep_design_bar_call <- function(exprs, group) {
    rhs <- as.call(list(as.name("+"), 0, .sep_add_calls(exprs)))
    as.call(list(as.name("|"),
                 rhs,
                 group))
}

.sep_spec_df <- function(x, what = "margin") {
    if (is.null(x)) return(NULL)
    if (is.matrix(x) || is.data.frame(x)) {
        ans <- as.data.frame(x, stringsAsFactors = FALSE)
    } else {
        ans <- data.frame(struc = unname(names(x)),
                          var = unname(x),
                          stringsAsFactors = FALSE)
    }
    if (!all(ans$struc %in% names(.sep_margin_registry))) {
        bad <- unique(ans$struc[!ans$struc %in% names(.sep_margin_registry)])
        stop("Unsupported separable() ", what, ": ", paste(bad, collapse = ", "))
    }
    ans
}

.sep_parse_spec <- function(x) {
    if (is.list(x) && !is.null(x$grid) && !is.null(x$margins)) {
        x$margins <- as.data.frame(x$margins, stringsAsFactors = FALSE)
        if (is.null(x$scale)) {
            x$scale <- list(mode = "auto", margins = NULL)
        } else if (!is.list(x$scale) || is.null(x$scale$mode)) {
            x$scale <- list(
                mode = "margin",
                margins = as.data.frame(x$scale, stringsAsFactors = FALSE)
            )
        }
        return(x)
    }
    stop("Internal separable() spec is missing or malformed.")
}

.sep_resolve_spec_id <- function(id, sepSpecs) {
    if (!is.numeric(id) || length(id) != 1L || is.na(id) ||
        is.null(sepSpecs)) {
        stop("Internal separable() spec id is missing or out of range.")
    }
    id <- as.integer(id)
    if (id < 1L || id > length(sepSpecs)) {
        stop("Internal separable() spec id is missing or out of range.")
    }
    id
}

.sep_spec_from_id_or_value <- function(x, sepSpecs) {
    if (is.numeric(x) && length(x) == 1L) {
        return(sepSpecs[[.sep_resolve_spec_id(x, sepSpecs)]])
    }
    x
}

.sep_margin_entry <- function(code, density_kind, n_scale, n_corr) {
    n_scale_fun <- if (is.function(n_scale)) n_scale else function(n) n_scale
    n_corr_fun <- if (is.function(n_corr)) n_corr else function(n) n_corr
    list(
        code = code,
        density_kind = density_kind,
        can_scale = is.function(n_scale) || as.integer(n_scale) > 0L,
        n_scale = function(n) as.integer(n_scale_fun(n)),
        n_corr = function(n) as.integer(n_corr_fun(n))
    )
}

.sep_margin_registry <- list(
    diag = .sep_margin_entry("diag", "diag", function(n) n, 0L),
    homdiag = .sep_margin_entry("homdiag", "diag", 1L, 0L),
    cs = .sep_margin_entry("cs", "dense_corr", function(n) n, 1L),
    homcs = .sep_margin_entry("homcs", "dense_corr", 1L, 1L),
    us = .sep_margin_entry("us", "dense_corr", function(n) n,
                           function(n) n * (n - 1L) / 2L),
    ar1 = .sep_margin_entry("ar1", "ar1", 0L, 1L),
    hetar1 = .sep_margin_entry("hetar1", "ar1", function(n) n, 1L),
    ou = .sep_margin_entry("ou", "spatial", 0L, 1L),
    exp = .sep_margin_entry("exp", "spatial", 0L, 1L),
    gau = .sep_margin_entry("gau", "spatial", 0L, 1L),
    mat = .sep_margin_entry("mat", "spatial", 0L, 2L),
    toep = .sep_margin_entry("toep", "toep", function(n) n,
                             function(n) n - 1L),
    homtoep = .sep_margin_entry("homtoep", "toep", 1L,
                                function(n) n - 1L)
)

.sep_density_kind_code <- c(
    dense_corr = 1L,
    ar1 = 2L,
    diag = 3L,
    spatial = 4L,
    toep = 5L
)

.sep_dispatch_code <- c(corr_corr = 1L)

.sep_scale_mode_code <- c(
    margin = 1L,          # one margin supplies absolute SDs
    global = 2L,          # one global scale, all margins correlation-only
    product = 3L,         # all eligible margin scales multiply
    selected_product = 4L    # selected margin scales multiply
)

.sep_corr_matrix_codes <- list(
    dense_corr = c("cs", "homcs", "us"),
    ar1 = c("ar1", "hetar1"),
    diag = c("diag", "homdiag"),
    spatial = c("ou", "exp", "gau", "mat"),
    toep = c("toep", "homtoep")
)

.sep_dispatch <- function(regs) {
    kinds <- vapply(regs, `[[`, character(1), "density_kind")
    codes <- vapply(regs, `[[`, character(1), "code")
    for (i in seq_along(regs)) {
        if (!kinds[[i]] %in% names(.sep_corr_matrix_codes) ||
            !codes[[i]] %in% .sep_corr_matrix_codes[[kinds[[i]]]]) {
            return(NA_character_)
        }
    }
    "corr_corr"
}

.sep_margin_label <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    paste0(x$struc, "(", x$var, ")", collapse = " x ")
}

.sep_scale_label <- function(scale) {
    if (is.null(scale)) return("NULL")
    if (is.null(scale$margins) || nrow(scale$margins) == 0L) {
        return(paste0(scale$mode, "()"))
    }
    paste0(scale$mode, "(",
           paste0(scale$margins$struc, "(", scale$margins$var, ")",
                  collapse = ", "),
           ")")
}

.sep_scale_info <- function(margins, regs, scale = NULL) {
    ## Resolve how absolute SD parameters enter the separable covariance.
    can_scale <- vapply(regs, `[[`, logical(1), "can_scale")
    scale_candidates <- which(can_scale)
    scale_mode <- if (is.null(scale)) "auto" else scale$mode

    if (identical(scale_mode, "auto")) {
        if (length(scale_candidates) == 0L) {
            stop("separable() margins ", .sep_margin_label(margins),
                 " define only a correlation product. Use scale = global() ",
                 "to add an overall scale.")
        }
        if (length(scale_candidates) > 1L) {
            stop("More than one separable() margin can carry scale in ",
                 .sep_margin_label(margins), ". Please specify the scale mode ",
                 "explicitly, for example scale = global(), scale = product(), ",
                 "or scale = ", margins$struc[scale_candidates[1]],
                 "(", margins$var[scale_candidates[1]], ").")
        }
        scale_mode <- "margin"
        scale_margin <- scale_candidates
    } else if (identical(scale_mode, "margin")) {
        scale_spec <- scale$margins
        if (nrow(scale_spec) != 1L) {
            stop("separable() scale must be a single margin call such as ",
                 "scale = us(0 + member).")
        }
        scale_margin <- which(margins$struc == scale_spec$struc &
                              margins$var == scale_spec$var)
        if (length(scale_margin) != 1L) {
            stop("separable() scale must match one of the specified margins, ",
                 "for example scale = us(0 + member) when us(0 + member) ",
                 "is a margin.")
        }
        if (!regs[[scale_margin]]$can_scale) {
            stop("separable() scale = ", scale_spec$struc, "(",
                 scale_spec$var, ") selects a correlation-only margin. ",
                 "Use a scale-capable margin or scale = global().")
        }
    } else if (identical(scale_mode, "global")) {
        scale_margin <- integer()
    } else if (identical(scale_mode, "product")) {
        scale_margin <- scale_candidates
        if (length(scale_margin) == 0L) {
            stop("separable() scale = product() needs at least one ",
                 "scale-capable margin.")
        }
    } else if (identical(scale_mode, "selected_product")) {
        scale_spec <- scale$margins
        margin_key <- paste(margins$struc, margins$var, sep = "\r")
        scale_key <- paste(scale_spec$struc, scale_spec$var, sep = "\r")
        scale_margin <- match(scale_key, margin_key)
        if (anyNA(scale_margin)) {
            i <- which(is.na(scale_margin))[[1]]
            stop("separable() scale margin ", scale_spec$struc[i], "(",
                 scale_spec$var[i], ") must match one of the specified ",
                 "margins.")
        }
        scale_ok <- vapply(regs[scale_margin], `[[`, logical(1), "can_scale")
        if (!all(scale_ok)) {
            i <- which(!scale_ok)[[1]]
            stop("separable() scale = ", scale_spec$struc[i], "(",
                 scale_spec$var[i], ") selects a correlation-only margin. ",
                 "Use a scale-capable margin or scale = global().")
        }
        if (anyDuplicated(scale_margin)) {
            stop("separable() scale margins must be unique.")
        }
    } else {
        stop("Unknown separable() scale mode: ", scale_mode)
    }

    list(
        mode = scale_mode,
        mode_code = as.integer(.sep_scale_mode_code[[scale_mode]]),
        spec = as.integer(scale_margin - 1L),
        margin = scale_margin
    )
}

.sep_stop_unsupported_dispatch <- function(margins, regs, scale = NULL) {
    ## Diagnose scale errors before reporting unsupported density combinations.
    .sep_scale_info(margins, regs, scale)
    stop("separable() frontend parsed ", .sep_margin_label(margins),
         ", but the backend currently only evaluates products among diag(), ",
         "homdiag(), ar1(), hetar1(), cs(), homcs(), us(), ou(), exp(), ",
         "gau(), mat(), toep(), and homtoep().")
}

.sep_spatial_info <- function(margins, spec, dims) {
    spatial <- margins$struc %in% c("ou", "exp", "gau", "mat")
    starts <- rep.int(-1L, length(dims))
    dists <- numeric()
    if (!any(spatial)) return(list(starts = starts, dists = dists))

    if (is.null(spec$margin_cnms)) {
        stop("separable() spatial margins require product-margin column names.")
    }
    for (i in which(spatial)) {
        coords <- tryCatch(suppressWarnings(parseNumLevels(spec$margin_cnms[[i]])),
                           error = function(e) {
                               stop("separable() spatial margins require ",
                                    "numeric coordinate levels, usually from ",
                                    "numFactor().", call. = FALSE)
                           })
        if (nrow(coords) != dims[[i]]) {
            stop("separable() spatial margin metadata does not match the ",
                 "product design.")
        }
        if (margins$struc[[i]] == "ou" && ncol(coords) != 1L) {
            stop("'ou' separable() margins are for 1D coordinates only.")
        }
        starts[[i]] <- length(dists)
        dists <- c(dists, as.vector(as.matrix(stats::dist(coords))))
    }
    list(starts = starts, dists = dists)
}

.sep_restruc_info <- function(spec, cnms, blksize) {
    ## R-side contract for currently supported separable terms.
    spec <- .sep_parse_spec(spec)

    if (!is.null(spec$dims)) {
        dims <- as.integer(spec$dims)
    } else {
        coords <- parseNumLevels(cnms)
        dims <- as.integer(apply(coords, 2, function(z) length(unique(z))))
    }
    if (prod(dims) != blksize)
        stop("separable() requires a complete rectangular product design.")

    margins <- .sep_spec_df(spec$margins)
    if (nrow(margins) != length(dims))
        stop("separable() margin metadata does not match the product design.")

    strucs <- margins$struc
    regs <- .sep_margin_registry[strucs]
    dispatch <- .sep_dispatch(regs)
    if (is.na(dispatch)) {
        .sep_stop_unsupported_dispatch(margins, regs, spec$scale)
    }
    if (!identical(margins$var, spec$grid)) {
        stop("The separable() margin variables must match the product design. ",
             "Use, for example, ",
             "separable(homcs(0 + member) %x% ar1(0 + time) | group).")
    }

    scale_info <- .sep_scale_info(margins, regs, spec$scale)

    scale_ntheta <- if (identical(scale_info$mode, "global")) {
        1L
    } else {
        sum(vapply(scale_info$margin,
                   function(i) regs[[i]]$n_scale(dims[[i]]),
                   integer(1)))
    }
    corr_ntheta <- sum(vapply(seq_along(regs), function(i) {
        regs[[i]]$n_corr(dims[[i]])
    }, integer(1)))
    ntheta <- scale_ntheta + corr_ntheta
    density_kind <- vapply(regs, `[[`, character(1), "density_kind")
    spatial_info <- .sep_spatial_info(margins, spec, dims)

    list(
        dims = dims,
        codes = as.integer(vapply(strucs, function(z) .valid_covstruct[[z]], numeric(1))),
        density_kinds = as.integer(.sep_density_kind_code[density_kind]),
        dispatch = as.integer(.sep_dispatch_code[dispatch]),
        scale_mode = scale_info$mode_code,
        scale_spec = scale_info$spec,
        dist_starts = spatial_info$starts,
        dists = spatial_info$dists,
        ntheta = as.integer(ntheta),
        density_kind = density_kind,
        margins = margins,
        scale = spec$scale
    )
}

.sep_flatten_product <- function(x) {
    if (is.call(x) && identical(.sep_call_name(x), "%x%") && length(x) == 3L) {
        return(c(.sep_flatten_product(x[[2]]), .sep_flatten_product(x[[3]])))
    }
    list(x)
}

.sep_parse_scale_arg <- function(scale) {
    if (is.null(scale)) return(list(mode = "auto", margins = NULL))

    if (is.call(scale)) {
        nm <- .sep_call_name(scale)
        if (identical(nm, "global") && length(scale) == 1L) {
            return(list(mode = nm, margins = NULL))
        }
        if (identical(nm, "product")) {
            args <- as.list(scale[-1])
            if (length(args) == 0L) {
                return(list(mode = "product", margins = NULL))
            }
            margins <- do.call(rbind, lapply(args, .sep_product_margin_spec))
            margins <- .sep_spec_df(margins, "scale margin")
            rownames(margins) <- NULL
            return(list(mode = "selected_product", margins = margins))
        }
    }

    list(mode = "margin",
         margins = .sep_spec_df(.sep_product_margin_spec(scale), "scale margin"))
}

.sep_make_product_spec <- function(bar_expr, scale = NULL) {
    ## Compile the public product syntax
    ##
    ##   separable(us(0 + role) %x% ar1(0 + day) | group, scale = us(0 + role))
    ##
    ## into an internal product-design representation.  Margins may contain
    ## multiple no-intercept columns; the current backend supports products of
    ## correlation-matrix margins.
    if (!is.call(bar_expr) || !identical(.sep_call_name(bar_expr), "|") ||
        length(bar_expr) != 3L) {
        stop("separable() product syntax must look like ",
             "separable(us(0 + role) %x% ar1(0 + day) | group, ...).")
    }
    margin_calls <- .sep_flatten_product(bar_expr[[2]])
    if (length(margin_calls) < 2L) {
        stop("separable() product syntax requires margins joined by %x%, ",
             "e.g. us(0 + role) %x% ar1(0 + day).")
    }
    margins <- do.call(rbind, lapply(margin_calls, .sep_product_margin_spec))
    rownames(margins) <- NULL
    if (anyDuplicated(margins$var)) {
        stop("separable() product margins must use distinct variables.")
    }
    if (!all(margins$struc %in% names(.sep_margin_registry))) {
        bad <- unique(margins$struc[!margins$struc %in% names(.sep_margin_registry)])
        stop("Unsupported separable() margin: ", paste(bad, collapse = ", "))
    }

    scale_spec <- .sep_parse_scale_arg(scale)

    structure(
        list(grid = margins$var,
             margins = margins,
             group = bar_expr[[3]],
             scale = scale_spec,
             grid_expr = .sep_design_bar_call(margins$expr, bar_expr[[3]])),
        class = "glmmTMB_separable_spec"
    )
}

.sep_margin_formula <- function(expr, env) {
    stats::as.formula(as.call(list(as.name("~"), expr)), env = env)
}

.sep_product_colnames <- function(cnms) {
    Reduce(function(a, b) as.vector(outer(a, b, paste, sep = ":")), cnms)
}

.sep_margin_matrix <- function(expr, fr, env) {
    f <- .sep_margin_formula(expr, env)
    X <- Matrix::sparse.model.matrix(f, data = fr)
    if ("(Intercept)" %in% colnames(X)) {
        stop("separable() product margins must be no-intercept formulas, ",
             "for example us(0 + member) %x% ar1(0 + time).")
    }
    X
}

.sep_sparse_rows <- function(X, n) {
    X <- methods::as(X, "TsparseMatrix")
    if (!length(X@x)) {
        return(rep(list(list(j = integer(0), x = numeric(0))), n))
    }
    rows <- split(data.frame(j = X@j + 1L, x = X@x), X@i + 1L)
    lapply(seq_len(n), function(i) {
        r <- rows[[as.character(i)]]
        if (is.null(r)) list(j = integer(0), x = numeric(0))
        else list(j = as.integer(r$j), x = as.numeric(r$x))
    })
}

.sep_row_kron_entries <- function(entries, dims) {
    ## First margin is fastest, matching expand.grid(), sepgrid(), and the C++
    ## array order used by the separable likelihood.
    ans <- list(j = 1L, x = 1)
    stride <- 1L
    for (m in seq_along(entries)) {
        e <- entries[[m]]
        if (!length(e$j) || !length(ans$j)) {
            return(list(j = integer(0), x = numeric(0)))
        }
        ans <- list(
            j = as.integer(as.vector(outer(ans$j, stride * (e$j - 1L), "+"))),
            x = as.vector(outer(ans$x, e$x, "*"))
        )
        stride <- stride * dims[[m]]
    }
    ans
}

.sep_build_product_reterm <- function(spec, fr, group, env) {
    margins <- .sep_spec_df(spec$margins)
    Xlist <- lapply(margins$expr, .sep_margin_matrix, fr = fr, env = env)
    dims <- vapply(Xlist, ncol, integer(1))
    cnms <- .sep_product_colnames(lapply(Xlist, colnames))

    g <- as.integer(group)
    if (anyNA(g)) stop("separable() grouping factor contains NA values.")
    n <- nrow(fr)
    p <- prod(dims)
    rows <- lapply(Xlist, .sep_sparse_rows, n = n)
    row_entries <- lapply(seq_len(n), function(r) {
        .sep_row_kron_entries(lapply(rows, `[[`, r), dims)
    })
    nnz <- lengths(lapply(row_entries, `[[`, "j"))
    obs <- rep.int(seq_len(n), nnz)
    jj <- unlist(lapply(row_entries, `[[`, "j"), use.names = FALSE)
    xx <- unlist(lapply(row_entries, `[[`, "x"), use.names = FALSE)
    Zt <- Matrix::sparseMatrix(
        i = (g[obs] - 1L) * p + jj,
        j = obs,
        x = xx,
        dims = c(nlevels(group) * p, n)
    )

    spec$dims <- dims
    spec$margin_cnms <- lapply(Xlist, colnames)
    spec$cnms <- cnms
    list(Zt = Zt, cnms = cnms, spec = spec)
}

.sep_replace_product_reterms <- function(reTrms, ss, sepSpecs, fr, env) {
    if (is.null(sepSpecs) || !length(sepSpecs)) return(list(reTrms = reTrms,
                                                            sepSpecs = sepSpecs))
    sep_pos <- which(ss$reTrmClasses == "separable")
    if (!length(sep_pos)) return(list(reTrms = reTrms, sepSpecs = sepSpecs))

    ## Run after smooth augmentation: by this point Ztlist positions match
    ## splitForm() term order, so separable terms can be replaced in place.
    assign <- attr(reTrms$flist, "assign")
    for (i in sep_pos) {
        id <- .sep_resolve_spec_id(eval(ss$reTrmAddArgs[[i]][[2]],
                                        envir = fr, enclos = env),
                                   sepSpecs)
        group <- reTrms$flist[[assign[i]]]
        repl <- .sep_build_product_reterm(sepSpecs[[id]], fr, group, env)
        reTrms$Ztlist[[i]] <- repl$Zt
        reTrms$cnms[[i]] <- repl$cnms
        sepSpecs[[id]] <- repl$spec
    }
    reTrms$Zt <- do.call(rbind, reTrms$Ztlist)
    reTrms$Gp <- cumsum(c(0L, vapply(reTrms$Ztlist, nrow, integer(1))))
    list(reTrms = reTrms, sepSpecs = sepSpecs)
}

.sep_reXterms <- function(spec, env) {
    margins <- .sep_spec_df(spec$margins)
    terms <- lapply(margins$expr, function(expr) {
        stats::terms(.sep_margin_formula(expr, env))
    })
    structure(list(margins = margins[, c("struc", "var"), drop = FALSE],
                   terms = terms,
                   cnms = spec$margin_cnms,
                   product_cnms = spec$cnms),
              class = "separable_reXterms")
}

.sep_replace_reXterms <- function(reXterms, ss, aa, sepSpecs, env) {
    if (is.null(sepSpecs) || !length(sepSpecs)) return(reXterms)
    sep_pos <- which(ss == "separable")
    for (i in sep_pos) {
        id <- .sep_resolve_spec_id(aa[[i]], sepSpecs)
        reXterms[[i]] <- .sep_reXterms(sepSpecs[[id]], env)
    }
    reXterms
}

.rewrite_separable_expr <- function(x, specs) {
    ## Walk the formula call tree and rewrite rich separable calls into the
    ## lower-level form that `reformulas::splitForm()` already understands:
    ##
    ##   separable(grid + 0 | group, <spec id>)
    ##
    ## The public syntax handled here is the product form:
    ##
    ##   separable(us(0 + member) %x% ar1(0 + time) | group,
    ##             scale = us(0 + member))
    ##
    ## If a formula already contains the lower-level list form, leave it alone.
    if (!is.call(x)) return(list(expr = x, specs = specs))
    if (identical(.sep_call_name(x), "separable")) {
        args <- as.list(x[-1])
        nms <- names(args)
        if (is.null(nms)) nms <- rep("", length(args))
        scale_i <- which(nms == "scale")
        if (length(scale_i) > 1L)
            stop("separable() accepts at most one scale argument.")
        scale_arg <- if (length(scale_i)) args[[scale_i]] else NULL
        named_i <- which(nzchar(nms) & nms != "scale")
        if (length(named_i) > 0L) {
            stop("separable() only accepts named argument scale.")
        }

        unnamed_args <- args[!nzchar(nms)]
        if (length(unnamed_args) == 1L) {
            spec <- .sep_make_product_spec(unnamed_args[[1]], scale = scale_arg)
            id <- length(specs) + 1L
            specs[[id]] <- spec
            return(list(expr = as.call(list(as.name("separable"),
                                            spec$grid_expr,
                                            as.integer(id))),
                        specs = specs))
        }

        ## Leave the older internal low-level form `separable(grid, list(...))`
        ## alone. New rewritten formulas use a numeric spec id instead.
        if (length(unnamed_args) == 2L && is.call(unnamed_args[[2]]) &&
            identical(.sep_call_name(unnamed_args[[2]]), "list") &&
            is.null(scale_arg)) {
            return(list(expr = x, specs = specs))
        }

        stop("separable() requires separable(",
             "margin1(0 + variable) %x% margin2(0 + variable) | group).")
    }
    for (i in seq_along(x)[-1]) {
        y <- .rewrite_separable_expr(x[[i]], specs)
        x[[i]] <- y$expr
        specs <- y$specs
    }
    list(expr = x, specs = specs)
}

rewrite_separable_formula <- function(f) {
    ## Apply the rewrite only to the RHS.  We call this early in `glmmTMB()` for
    ## conditional, zero-inflation, and dispersion formulas so all downstream
    ## machinery sees a normal random-effect special.
    if (!inherits(f, "formula")) return(f)
    y <- .rewrite_separable_expr(f[[length(f)]], list())
    f[[length(f)]] <- y$expr
    if (length(y$specs) > 0L) {
        attr(f, "separable_specs") <- y$specs
    }
    f
}
