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

## The helpers below parse the public product syntax
##
##   separable(us(0 + member) %x% ar1(0 + time) | group,
##             scale = us(0 + member))
##
## into an internal spec.  Separable terms are split with `splitForm()` like
## other covariance specials, but are excluded from `mkReTrms()` because the
## marginal covariance calls are not ordinary model-matrix expressions.
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

.sep_formula_specs <- function(f) {
    if (!inherits(f, "formula")) return(list())
    ss <- reformulas::splitForm(f, specials = c(names(.valid_covstruct), "s"))
    .sep_specs_from_split(ss)
}

.sep_margin_extra_varnames <- function(margins) {
    unlist(Map(function(struc, extra) {
        frame_args <- .sep_margin_registry[[struc]]$extra$frame_args
        unlist(lapply(extra[frame_args], all.vars), use.names = FALSE)
    }, margins$struc, margins$extra), use.names = FALSE)
}

.sep_margin_varnames <- function(..., include_group = FALSE) {
    forms <- list(...)
    vars <- unlist(lapply(forms, function(f) {
        specs <- .sep_formula_specs(f)
        unlist(lapply(specs, function(spec) {
            ans <- unlist(lapply(spec$margins$expr, all.vars), use.names = FALSE)
            ans <- c(ans, .sep_margin_extra_varnames(spec$margins))
            if (include_group) ans <- c(ans, all.vars(spec$group))
            unique(ans)
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
    if (!is.call(x) || length(x) < 2L)
        stop("separable() product margins must look like us(0 + role) ",
             "or ar1(0 + day).")
    args <- as.list(x[-1])
    data.frame(struc = .sep_deparse(x[[1]]),
               var = .sep_product_margin_label(args[[1]]),
               expr = I(list(args[[1]])),
               extra = I(list(args[-1])),
               stringsAsFactors = FALSE)
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
    if (is.null(ans$extra)) {
        ans$extra <- I(rep(list(list()), nrow(ans)))
    }
    if (!all(ans$struc %in% names(.sep_margin_registry))) {
        bad <- unique(ans$struc[!ans$struc %in% names(.sep_margin_registry)])
        stop("Unsupported separable() ", what, ": ", paste(bad, collapse = ", "))
    }
    ans
}

.sep_validate_margin_extras <- function(margins, what = "margin") {
    margins <- .sep_spec_df(margins, what)
    n_extra <- lengths(margins$extra)
    n_expected <- vapply(margins$struc, function(z) {
        .sep_margin_registry[[z]]$extra$n
    }, integer(1))
    bad <- which(n_extra != n_expected)
    if (length(bad)) {
        i <- bad[[1]]
        stop("separable() ", what, " ", margins$struc[i], "(",
             margins$var[i], ") takes ", n_expected[[i]],
             " extra argument", if (n_expected[[i]] == 1L) "" else "s",
             ", but got ", n_extra[[i]], ".")
    }
    margins
}

.sep_parse_spec <- function(x) {
    if (is.list(x) && !is.null(x$grid) && !is.null(x$margins)) {
        x$margins <- .sep_spec_df(x$margins)
        if (is.null(x$scale)) {
            x$scale <- list(mode = "auto", margins = NULL)
        } else if (!is.list(x$scale) || is.null(x$scale$mode)) {
            x$scale <- list(
                mode = "margin",
                margins = .sep_spec_df(x$scale)
            )
        }
        return(x)
    }
    stop("Internal separable() spec is missing or malformed.")
}

.sep_nfun <- function(n) {
    if (is.function(n)) return(function(d) as.integer(n(d)))
    force(n)
    function(d) as.integer(n)
}

.sep_theta <- function(...) {
    blocks <- list(...)
    if (length(blocks) == 1L && is.null(names(blocks))) {
        names(blocks) <- "corr"
    }
    list(blocks = Map(function(name, n) {
        list(name = name, n = .sep_nfun(n))
    }, names(blocks), blocks))
}

.sep_scale <- function(kind = c("none", "homogeneous", "heterogeneous"),
                       n = 0L, auto = !identical(kind, "none"),
                       fixed = FALSE) {
    kind <- match.arg(kind)
    list(kind = kind,
         can_scale = !identical(kind, "none"),
         can_auto_scale = isTRUE(auto),
         fixed_scale = isTRUE(fixed),
         n = .sep_nfun(n))
}

.sep_metadata <- function(dist_coord_dim = NULL) {
    needs_dist <- !is.null(dist_coord_dim)
    list(needs_dist = needs_dist,
         dist_coord_dim = if (needs_dist) as.integer(dist_coord_dim) else NA_integer_)
}

.sep_payload <- function(kind, arg, frame = FALSE) {
    list(kind = kind, arg = as.integer(arg), frame = isTRUE(frame))
}

.sep_extra <- function(n = 0L, frame_args = integer(), payloads = list(),
                       cov_args = integer()) {
    n <- as.integer(n)
    frame_args <- as.integer(frame_args)
    if (length(cov_args)) {
        payloads <- c(payloads, lapply(as.integer(cov_args), function(i) {
            .sep_payload("cov_matrix", i)
        }))
    }
    if (length(payloads) && !is.list(payloads[[1]])) {
        payloads <- list(payloads)
    }
    payloads <- lapply(payloads, function(x) {
        .sep_payload(x$kind, x$arg, x$frame)
    })
    payload_args <- vapply(payloads, `[[`, integer(1), "arg")
    frame_args <- sort(unique(c(frame_args,
                                payload_args[vapply(payloads, `[[`,
                                                     logical(1), "frame")])))
    list(n = n, frame_args = frame_args, payloads = payloads)
}

.sep_margin_entry <- function(scale, theta = .sep_theta(),
                              metadata = .sep_metadata(),
                              extra = .sep_extra(),
                              validate_payload = NULL,
                              require_grid_columns = FALSE) {
    list(
        scale = scale,
        theta = theta,
        metadata = metadata,
        extra = extra,
        validate_payload = validate_payload,
        require_grid_columns = require_grid_columns
    )
}

.sep_cov_matrix_extra <- .sep_extra(n = 1L, cov_args = 1L)

.sep_validate_propto_payload <- function(payload, value, cnms, expr, env) {
    if (!identical(payload$kind, "cov_matrix")) return(value)
    if (!is.numeric(value)) {
        stop("separable() propto() margin expects a numeric matrix ",
             "extra argument.", call. = FALSE)
    }
    tryCatch(
        checkProptoNames(aa = value, cnms = cnms,
                         reXtrm = .sep_margin_formula(expr, env)),
        error = function(e) {
            stop("separable() propto() margin ", conditionMessage(e),
                 call. = FALSE)
        })
    value
}

.sep_margin_registry <- local({
    n_dim <- function(n) n
    n_pairs <- function(n) n * (n - 1L) / 2L
    n_lags <- function(n) n - 1L
    het <- .sep_scale("heterogeneous", n_dim)
    hom <- .sep_scale("homogeneous", 1L)
    spatial <- .sep_metadata(NA)

    ans <- list(
        diag = .sep_margin_entry(het),
        homdiag = .sep_margin_entry(hom),
        cs = .sep_margin_entry(het, .sep_theta(corr = 1L)),
        homcs = .sep_margin_entry(hom, .sep_theta(corr = 1L)),
        us = .sep_margin_entry(het, .sep_theta(corr = n_pairs)),
        ar1 = .sep_margin_entry(.sep_scale("homogeneous", 1L, auto = FALSE),
                                .sep_theta(corr = 1L),
                                require_grid_columns = TRUE),
        hetar1 = .sep_margin_entry(het, .sep_theta(corr = 1L),
                                   require_grid_columns = TRUE),
        ou = .sep_margin_entry(hom, .sep_theta(decay = 1L),
                               .sep_metadata(1L),
                               require_grid_columns = TRUE),
        exp = .sep_margin_entry(hom, .sep_theta(range = 1L), spatial,
                                require_grid_columns = TRUE),
        gau = .sep_margin_entry(hom, .sep_theta(range = 1L), spatial,
                                require_grid_columns = TRUE),
        mat = .sep_margin_entry(hom,
                                .sep_theta(range = 1L, smoothness = 1L),
                                spatial,
                                require_grid_columns = TRUE),
        toep = .sep_margin_entry(het, .sep_theta(corr = n_lags),
                                 require_grid_columns = TRUE),
        homtoep = .sep_margin_entry(hom, .sep_theta(corr = n_lags),
                                    require_grid_columns = TRUE),
        propto = .sep_margin_entry(hom,
                                   extra = .sep_cov_matrix_extra,
                                   validate_payload =
                                       .sep_validate_propto_payload),
        equalto = .sep_margin_entry(.sep_scale("none", fixed = TRUE),
                                    extra = .sep_cov_matrix_extra)
    )
    for (code in names(ans)) ans[[code]]$code <- code
    ans
})

.sep_scale_kind_code <- c(
    none = 0L,
    homogeneous = 1L,
    heterogeneous = 2L
)

.sep_theta_block_kind_code <- c(
    global_scale = 1L,
    scale = 2L,
    corr = 3L,
    range = 4L,
    smoothness = 5L,
    decay = 6L
)

.sep_matrix_payload_kind_code <- c(
    distance = 1L,
    cov_matrix = 2L
)

.sep_margin_label <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    paste0(x$struc, "(", x$var, ")", collapse = " x ")
}

.sep_margin_key <- function(x) {
    x <- .sep_spec_df(x)
    extra <- vapply(x$extra, function(z) {
        paste(vapply(z, .sep_deparse, character(1)), collapse = "\r")
    }, character(1))
    paste(x$struc, x$var, extra, sep = "\r")
}

.sep_match_scale_margins <- function(scale_spec, margins, regs, single = FALSE,
                                     unique = FALSE) {
    if (single && nrow(scale_spec) != 1L)
        stop("separable() scale must be a single margin call such as ",
             "scale = us(0 + member).")
    margin_key <- .sep_margin_key(margins)
    scale_key <- .sep_margin_key(scale_spec)
    scale_margin <- if (single) which(margin_key == scale_key) else
        match(scale_key, margin_key)
    if (single && length(scale_margin) != 1L)
        stop("separable() scale must match one of the specified margins, ",
             "for example scale = us(0 + member) when us(0 + member) ",
             "is a margin.")
    if (anyNA(scale_margin)) {
        i <- which(is.na(scale_margin))[[1]]
        stop("separable() scale margin ", scale_spec$struc[i], "(",
             scale_spec$var[i], ") must match one of the specified ",
             "margins.")
    }
    bad_scale <- which(!vapply(regs[scale_margin],
                               function(x) x$scale$can_scale, logical(1)))
    if (length(bad_scale)) {
        i <- bad_scale[[1]]
        stop("separable() scale = ", scale_spec$struc[i], "(",
             scale_spec$var[i], ") selects a correlation-only margin. ",
             "Use a scale-capable margin or scale = global().")
    }
    if (unique && anyDuplicated(scale_margin))
        stop("separable() scale margins must be unique.")
    scale_margin
}

.sep_scale_info <- function(margins, regs, scale = NULL) {
    ## Resolve how absolute SD parameters enter the separable covariance.
    can_scale <- vapply(regs, function(x) x$scale$can_scale, logical(1))
    can_auto_scale <- vapply(regs, function(x) x$scale$can_auto_scale,
                              logical(1))
    fixed_scale <- vapply(regs, function(x) x$scale$fixed_scale, logical(1))
    scale_candidates <- which(can_auto_scale)
    scale_mode <- if (is.null(scale)) "auto" else scale$mode

    if (identical(scale_mode, "auto")) {
        if (length(scale_candidates) == 0L) {
            if (any(fixed_scale)) {
                scale_mode <- "none"
                scale_margin <- integer()
            } else {
                stop("separable() margins ", .sep_margin_label(margins),
                     " have no unambiguous scale margin. Use scale = global() ",
                     "or choose a scale margin explicitly.")
            }
        } else if (length(scale_candidates) > 1L) {
            stop("More than one separable() margin can carry scale in ",
                 .sep_margin_label(margins), ". Please specify the scale mode ",
                 "explicitly, for example scale = global(), scale = product(), ",
                 "or scale = ", margins$struc[scale_candidates[1]],
                 "(", margins$var[scale_candidates[1]], ").")
        } else {
            scale_mode <- "margin"
            scale_margin <- scale_candidates
        }
    } else if (identical(scale_mode, "margin")) {
        scale_margin <- .sep_match_scale_margins(scale$margins, margins, regs,
                                                 single = TRUE)
    } else if (identical(scale_mode, "global")) {
        scale_margin <- integer()
    } else if (identical(scale_mode, "product")) {
        scale_margin <- which(can_scale)
        if (length(scale_margin) == 0L) {
            stop("separable() scale = product() needs at least one ",
                 "scale-capable margin.")
        }
    } else if (identical(scale_mode, "selected_product")) {
        scale_margin <- .sep_match_scale_margins(scale$margins, margins, regs,
                                                 unique = TRUE)
    } else if (identical(scale_mode, "none")) {
        scale_margin <- integer()
    } else {
        stop("Unknown separable() scale mode: ", scale_mode)
    }

    list(
        mode = scale_mode,
        spec = as.integer(scale_margin - 1L),
        margin = scale_margin
    )
}

.sep_registry_matrix_payloads <- function(regs) {
    lapply(regs, function(reg) {
        payloads <- reg$extra$payloads
        keep <- vapply(payloads, function(x) {
            x$kind %in% names(.sep_matrix_payload_kind_code)
        }, logical(1))
        payloads[keep]
    })
}

.sep_matrix_payload_info <- function(regs, spec, dims) {
    needs_dist <- vapply(regs, function(x) x$metadata$needs_dist, logical(1))
    matrix_payloads <- .sep_registry_matrix_payloads(regs)
    ## This full matrix-payload kind table (sepMatrixPayloadKinds) indexes
    ## per-margin starts, not just kinds present in the current separable term.
    n_kind <- length(.sep_matrix_payload_kind_code)
    starts <- rep.int(-1L, length(dims) * n_kind)
    values <- numeric()

    add_matrix <- function(i, kind, mat) {
        if (nrow(mat) != dims[[i]] || ncol(mat) != dims[[i]]) {
            stop("separable() ", kind, " matrix metadata does not match ",
                 "the product design.")
        }
        k <- match(kind, names(.sep_matrix_payload_kind_code))
        starts[(i - 1L) * n_kind + k] <<- length(values)
        values <<- c(values, as.vector(mat))
    }

    if (any(needs_dist) && is.null(spec$margin_cnms)) {
        stop("separable() spatial margins require product-margin column names.")
    }
    for (i in which(needs_dist)) {
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
        coord_dim <- regs[[i]]$metadata$dist_coord_dim
        if (!is.na(coord_dim) && ncol(coords) != coord_dim) {
            stop("'", regs[[i]]$code, "' separable() margins are for ",
                 coord_dim, "D coordinates only.")
        }
        add_matrix(i, "distance", as.matrix(stats::dist(coords)))
    }
    for (i in seq_along(matrix_payloads)) {
        for (payload in matrix_payloads[[i]]) {
            mat <- spec$margin_payloads[[i]][[payload$kind]]
            if (is.null(mat)) {
                stop("separable() ", payload$kind, " margin is missing ",
                     "matrix metadata.")
            }
            add_matrix(i, payload$kind, mat)
        }
    }

    list(kinds = as.integer(.sep_matrix_payload_kind_code),
         starts = as.integer(starts),
         values = values)
}

.sep_theta_block_layout <- function(regs, dims, scale_info) {
    rows <- list()
    pos <- 0L

    add_block <- function(margin, kind, n) {
        if (!kind %in% names(.sep_theta_block_kind_code)) {
            stop("Unsupported separable() theta block kind: ", kind)
        }
        rows[[length(rows) + 1L]] <<- data.frame(
            margin = as.integer(margin),
            kind = as.integer(.sep_theta_block_kind_code[[kind]]),
            start = as.integer(pos),
            length = as.integer(n),
            stringsAsFactors = FALSE
        )
        pos <<- pos + as.integer(n)
    }

    if (identical(scale_info$mode, "global")) {
        add_block(-1L, "global_scale", 1L)
    }
    for (i in seq_along(regs)) {
        if (i %in% scale_info$margin) {
            add_block(i - 1L, "scale", regs[[i]]$scale$n(dims[[i]]))
        }
        for (block in regs[[i]]$theta$blocks) {
            add_block(i - 1L, block$name, block$n(dims[[i]]))
        }
    }

    ans <- if (length(rows)) do.call(rbind, rows) else
        data.frame(margin = integer(), kind = integer(),
                   start = integer(), length = integer())
    as.matrix(ans)
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
    if (!identical(margins$var, spec$grid)) {
        stop("The separable() margin variables must match the product design. ",
             "Use, for example, ",
             "separable(homcs(0 + member) %x% ar1(0 + time) | group).")
    }

    scale_info <- .sep_scale_info(margins, regs, spec$scale)

    theta_layout <- .sep_theta_block_layout(regs, dims, scale_info)
    ntheta <- sum(theta_layout[, "length"])
    scale_kind <- vapply(regs, function(x) x$scale$kind, character(1))
    matrix_info <- .sep_matrix_payload_info(regs, spec, dims)

    list(
        tmb = list(
            sepDims = dims,
            sepMarginStruc = strucs,
            sepMarginVars = margins$var,
            sepCodes = as.integer(vapply(strucs, function(z) .valid_covstruct[[z]],
                                         numeric(1))),
            sepScaleKinds = as.integer(.sep_scale_kind_code[scale_kind]),
            sepScaleSpec = scale_info$spec,
            sepThetaBlocks = theta_layout,
            sepMatrixPayloadKinds = matrix_info$kinds,
            sepMatrixPayloadStarts = matrix_info$starts,
            sepMatrixPayloadValues = matrix_info$values
        ),
        ntheta = as.integer(ntheta)
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
            margins <- .sep_validate_margin_extras(margins, "scale margin")
            rownames(margins) <- NULL
            return(list(mode = "selected_product", margins = margins))
        }
    }

    list(mode = "margin",
         margins = .sep_validate_margin_extras(.sep_product_margin_spec(scale),
                                               "scale margin"))
}

.sep_make_product_spec <- function(bar_expr, scale = NULL) {
    ## Parse the public product syntax
    ##
    ##   separable(us(0 + role) %x% ar1(0 + day) | group, scale = us(0 + role))
    ##
    ## into an internal product-design representation.  Margins may contain
    ## multiple no-intercept columns; the current backend supports registered
    ## margins that can be represented as marginal SDs and correlations.
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
    margins <- .sep_validate_margin_extras(margins)
    rownames(margins) <- NULL
    if (anyDuplicated(margins$var)) {
        stop("separable() product margins must use distinct variables.")
    }

    scale_spec <- .sep_parse_scale_arg(scale)

    structure(
        list(grid = margins$var,
             margins = margins,
             group = bar_expr[[3]],
             scale = scale_spec),
        class = "glmmTMB_separable_spec"
    )
}

.sep_scale_arg_from_add_arg <- function(add_arg) {
    if (!is.call(add_arg) || !identical(.sep_call_name(add_arg), "separable")) {
        stop("Internal separable() add-argument is malformed.")
    }
    args <- as.list(add_arg[-1])
    nms <- names(args)
    if (is.null(nms)) nms <- rep("", length(args))
    scale_i <- which(nms == "scale")
    if (length(scale_i) > 1L)
        stop("separable() accepts at most one scale argument.")
    named_i <- which(nzchar(nms) & nms != "scale")
    if (length(named_i) > 0L) {
        stop("separable() only accepts named argument scale.")
    }
    unnamed_args <- args[!nzchar(nms)]
    if (length(unnamed_args) > 0L) {
        stop("separable() requires separable(",
             "margin1(0 + variable) %x% margin2(0 + variable) | group).")
    }
    if (length(scale_i)) args[[scale_i]] else NULL
}

.sep_make_product_spec_from_split <- function(bar_expr, add_arg) {
    .sep_make_product_spec(bar_expr, .sep_scale_arg_from_add_arg(add_arg))
}

.sep_specs_from_split <- function(ss) {
    sep_pos <- which(ss$reTrmClasses == "separable")
    specs <- vector("list", length(ss$reTrmClasses))
    for (i in sep_pos) {
        specs[[i]] <- .sep_make_product_spec_from_split(ss$reTrmFormulas[[i]],
                                                        ss$reTrmAddArgs[[i]])
    }
    specs[!vapply(specs, is.null, logical(1))]
}

.sep_margin_formula <- function(expr, env) {
    stats::as.formula(as.call(list(as.name("~"), expr)), env = env)
}

.sep_product_colnames <- function(cnms) {
    Reduce(function(a, b) as.vector(outer(a, b, paste, sep = ":")), cnms)
}

.sep_margin_matrix <- function(expr, fr, env) {
    f <- .sep_margin_formula(expr, env)
    ## Matrix::sparse.model.matrix() treats a terms-bearing data frame as a
    ## prebuilt model frame, which rejects transformed margin expressions.
    attr(fr, "terms") <- NULL
    X <- Matrix::sparse.model.matrix(f, data = fr)
    if ("(Intercept)" %in% colnames(X)) {
        stop("separable() product margins must be no-intercept formulas, ",
             "for example us(0 + member) %x% ar1(0 + time).")
    }
    X
}

.sep_eval_extra <- function(expr, fr, env) {
    tryCatch(eval(expr, envir = fr, enclos = env),
             error = function(e) {
                 stop("can't evaluate separable() margin argument ",
                      sQuote(.sep_deparse(expr)), call. = FALSE)
             })
}

.sep_check_cov_payload <- function(x, cnms, struc) {
    if (identical(struc, "equalto")) {
        tryCatch(checkEqualto(aa = x, cnms = cnms),
                 error = function(e) {
                     stop("separable() equalto() margin ",
                          conditionMessage(e), call. = FALSE)
                 })
        return(x)
    }
    if (identical(struc, "propto")) {
        return(x)
    }
    if (!is.matrix(x) || !is.numeric(x) ||
        nrow(x) != ncol(x) || nrow(x) != length(cnms)) {
        stop("separable() ", struc, "() margin expects a numeric square ",
             "matrix matching the margin columns.", call. = FALSE)
    }
    x
}

.sep_check_payload <- function(payload, value, cnms, struc) {
    switch(payload$kind,
           cov_matrix = .sep_check_cov_payload(value, cnms, struc),
           stop("Unsupported separable() payload kind: ", payload$kind,
                call. = FALSE))
}

.sep_resolve_margin_payloads <- function(margins, Xlist, fr, env) {
    payloads <- vector("list", nrow(margins))
    for (i in seq_len(nrow(margins))) {
        reg <- .sep_margin_registry[[margins$struc[i]]]
        specs <- reg$extra$payloads
        if (!length(specs)) next
        cnms <- colnames(Xlist[[i]])
        payloads[[i]] <- setNames(lapply(specs, function(payload) {
            value <- .sep_eval_extra(margins$extra[[i]][[payload$arg]], fr, env)
            value <- .sep_check_payload(payload, value, cnms, margins$struc[i])
            if (is.function(reg$validate_payload)) {
                value <- reg$validate_payload(payload, value, cnms,
                                              margins$expr[[i]], env)
            }
            value
        }), vapply(specs, `[[`, character(1), "kind"))
    }
    payloads
}

.sep_validate_margin_matrices <- function(margins, Xlist) {
    regs <- .sep_margin_registry[margins$struc]
    needs_grid_columns <- vapply(regs, `[[`, logical(1),
                                 "require_grid_columns")
    bad <- which(needs_grid_columns &
                 vapply(Xlist, ncol, integer(1)) == 1L)
    if (length(bad)) {
        i <- bad[[1]]
        stop("separable() ", margins$struc[i], "() margin must define ",
             "one model-matrix column per level/coordinate; use a factor ",
             "or numFactor() rather than a numeric covariate.",
             call. = FALSE)
    }
}

.sep_build_product_reterm <- function(spec, fr, group, env) {
    margins <- .sep_spec_df(spec$margins)
    Xlist <- lapply(margins$expr, .sep_margin_matrix, fr = fr, env = env)
    .sep_validate_margin_matrices(margins, Xlist)
    dims <- vapply(Xlist, ncol, integer(1))
    cnms <- .sep_product_colnames(lapply(Xlist, colnames))

    g <- as.integer(group)
    if (anyNA(g)) stop("separable() grouping factor contains NA values.")
    n <- nrow(fr)
    product_t <- Reduce(function(prod_t, X) {
        Matrix::KhatriRao(Matrix::t(X), prod_t)
    }, Xlist[-1L], init = Matrix::t(Xlist[[1L]]))
    groups_t <- Matrix::sparseMatrix(i = g, j = seq_len(n), x = 1,
                                     dims = c(nlevels(group), n))
    Zt <- Matrix::KhatriRao(groups_t, product_t)

    spec$dims <- dims
    spec$margin_cnms <- lapply(Xlist, colnames)
    spec$margin_payloads <- .sep_resolve_margin_payloads(margins, Xlist, fr, env)
    spec$cnms <- cnms
    list(Zt = Zt, cnms = cnms, spec = spec)
}

.sep_group_factor <- function(expr, fr, env) {
    g <- eval(expr, envir = fr, enclos = env)
    if (anyNA(g)) stop("separable() grouping factor contains NA values.")
    factor(g)
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
