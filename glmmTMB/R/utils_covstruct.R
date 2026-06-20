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

.sep_theta <- function(...) {
    blocks <- list(...)
    nms <- names(blocks)
    if (length(blocks) == 1L && (is.null(nms) || !nzchar(nms))) {
        nms <- "corr"
    }
    if (length(blocks) > 0L &&
        (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms))) {
        stop("Malformed separable() theta contract.")
    }
    blocks <- Map(function(name, n) {
        n_fun <- if (is.function(n)) n else function(d) n
        list(name = name, n = function(d) as.integer(n_fun(d)))
    }, nms, blocks)
    list(blocks = blocks)
}

.sep_scale <- function(kind = c("none", "homogeneous", "heterogeneous"),
                       n = 0L, auto = NULL, fixed = FALSE) {
    kind <- match.arg(kind)
    n_fun <- if (is.function(n)) n else function(d) n
    can_scale <- !identical(kind, "none")
    if (is.null(auto)) auto <- can_scale
    fixed <- isTRUE(fixed)
    list(kind = kind,
         can_scale = can_scale,
         can_auto_scale = auto,
         fixed_scale = fixed,
         n = function(d) as.integer(n_fun(d)))
}

.sep_metadata <- function(dist_coord_dim = NULL) {
    needs_dist <- !is.null(dist_coord_dim)
    if (needs_dist && !(length(dist_coord_dim) == 1L &&
                        (is.na(dist_coord_dim) || dist_coord_dim >= 1L))) {
        stop("Unsupported separable() distance metadata.")
    }
    list(needs_dist = needs_dist,
         dist_coord_dim = if (needs_dist) as.integer(dist_coord_dim) else NA_integer_)
}

.sep_payload <- function(kind, arg, frame = FALSE) {
    if (!is.character(kind) || length(kind) != 1L || !nzchar(kind) ||
        length(arg) != 1L || is.na(arg) || arg < 1L ||
        !is.logical(frame) || length(frame) != 1L || is.na(frame)) {
        stop("Malformed separable() payload contract.")
    }
    list(kind = kind, arg = as.integer(arg), frame = frame)
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
        if (is.null(x$kind) || is.null(x$arg) || is.null(x$frame)) {
            stop("Malformed separable() payload contract.")
        }
        .sep_payload(x$kind, x$arg, x$frame)
    })
    if (length(n) != 1L || n < 0L ||
        any(frame_args < 1L | frame_args > n)) {
        stop("Malformed separable() extra-argument contract.")
    }
    payload_args <- vapply(payloads, `[[`, integer(1), "arg")
    if (any(payload_args < 1L | payload_args > n)) {
        stop("Malformed separable() payload contract.")
    }
    frame_args <- sort(unique(c(frame_args,
                                payload_args[vapply(payloads, `[[`,
                                                     logical(1), "frame")])))
    list(n = n, frame_args = frame_args, payloads = payloads)
}

.sep_margin_entry <- function(code, builder, scale, theta,
                              metadata = .sep_metadata(),
                              extra = .sep_extra()) {
    if (is.null(scale$kind) || !is.function(scale$n) ||
        is.null(scale$can_scale) || is.null(scale$can_auto_scale) ||
        is.null(scale$fixed_scale)) {
        stop("Malformed separable() scale contract for ", code)
    }
    if (is.null(theta$blocks)) {
        stop("Malformed separable() theta contract for ", code)
    }
    if (is.null(metadata$needs_dist) || is.null(metadata$dist_coord_dim)) {
        stop("Malformed separable() metadata contract for ", code)
    }
    if (is.null(extra$n) || is.null(extra$frame_args) ||
        is.null(extra$payloads)) {
        stop("Malformed separable() extra-argument contract for ", code)
    }
    if (!scale$kind %in% c("none", "homogeneous", "heterogeneous")) {
        stop("Unknown separable() scale kind for ", code, ": ", scale$kind)
    }
    list(
        code = code,
        builder = builder,
        scale = scale,
        theta = theta,
        metadata = metadata,
        extra = extra
    )
}

.sep_cov_matrix_extra <- .sep_extra(
    n = 1L,
    payloads = .sep_payload("cov_matrix", 1L)
)

.sep_margin_registry <- list(
    diag = .sep_margin_entry("diag", "diag",
                             .sep_scale("heterogeneous", function(n) n),
                             .sep_theta()),
    homdiag = .sep_margin_entry("homdiag", "diag",
                                .sep_scale("homogeneous", 1L),
                                .sep_theta()),
    cs = .sep_margin_entry("cs", "dense_corr",
                           .sep_scale("heterogeneous", function(n) n),
                           .sep_theta(corr = 1L)),
    homcs = .sep_margin_entry("homcs", "dense_corr",
                              .sep_scale("homogeneous", 1L),
                              .sep_theta(corr = 1L)),
    us = .sep_margin_entry("us", "dense_corr",
                           .sep_scale("heterogeneous", function(n) n),
                           .sep_theta(corr = function(n) n * (n - 1L) / 2L)),
    ar1 = .sep_margin_entry("ar1", "ar1",
                            .sep_scale("homogeneous", 1L, auto = FALSE),
                            .sep_theta(corr = 1L)),
    hetar1 = .sep_margin_entry("hetar1", "ar1",
                               .sep_scale("heterogeneous", function(n) n),
                               .sep_theta(corr = 1L)),
    ou = .sep_margin_entry("ou", "spatial",
                           .sep_scale("homogeneous", 1L),
                           .sep_theta(decay = 1L),
                           .sep_metadata(dist_coord_dim = 1L)),
    exp = .sep_margin_entry("exp", "spatial",
                            .sep_scale("homogeneous", 1L),
                            .sep_theta(range = 1L),
                            .sep_metadata(dist_coord_dim = NA)),
    gau = .sep_margin_entry("gau", "spatial",
                            .sep_scale("homogeneous", 1L),
                            .sep_theta(range = 1L),
                            .sep_metadata(dist_coord_dim = NA)),
    mat = .sep_margin_entry("mat", "spatial",
                            .sep_scale("homogeneous", 1L),
                            .sep_theta(range = 1L, smoothness = 1L),
                            .sep_metadata(dist_coord_dim = NA)),
    toep = .sep_margin_entry("toep", "toep",
                             .sep_scale("heterogeneous", function(n) n),
                             .sep_theta(corr = function(n) n - 1L)),
    homtoep = .sep_margin_entry("homtoep", "toep",
                                .sep_scale("homogeneous", 1L),
                                .sep_theta(corr = function(n) n - 1L)),
    propto = .sep_margin_entry("propto", "fixed_cov",
                               .sep_scale("homogeneous", 1L),
                               .sep_theta(),
                               extra = .sep_cov_matrix_extra),
    equalto = .sep_margin_entry("equalto", "fixed_cov",
                                .sep_scale("none", fixed = TRUE),
                                .sep_theta(),
                                extra = .sep_cov_matrix_extra)
)

## Each supported separable builder kind must have a C++ margin builder that
## returns a standard-deviation vector and a correlation matrix.
.sep_builder_kind_code <- c(
    dense_corr = 1L,
    ar1 = 2L,
    diag = 3L,
    spatial = 4L,
    toep = 5L,
    fixed_cov = 6L
)

.sep_dispatch_code <- c(corr_matrix_product = 1L)

.sep_scale_mode_code <- c(
    none = 0L,            # fixed scales supplied by one or more margins
    margin = 1L,          # one margin supplies absolute SDs
    global = 2L,          # one global scale, all margins correlation-only
    product = 3L,         # all eligible margin scales multiply
    selected_product = 4L    # selected margin scales multiply
)

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

.sep_dispatch <- function(regs) {
    kinds <- vapply(regs, `[[`, character(1), "builder")
    if (!all(kinds %in% names(.sep_builder_kind_code))) return(NA_character_)
    "corr_matrix_product"
}

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
        scale_spec <- scale$margins
        if (nrow(scale_spec) != 1L) {
            stop("separable() scale must be a single margin call such as ",
                 "scale = us(0 + member).")
        }
        scale_margin <- which(.sep_margin_key(margins) ==
                              .sep_margin_key(scale_spec))
        if (length(scale_margin) != 1L) {
            stop("separable() scale must match one of the specified margins, ",
                 "for example scale = us(0 + member) when us(0 + member) ",
                 "is a margin.")
        }
        if (!regs[[scale_margin]]$scale$can_scale) {
            stop("separable() scale = ", scale_spec$struc, "(",
                 scale_spec$var, ") selects a correlation-only margin. ",
                 "Use a scale-capable margin or scale = global().")
        }
    } else if (identical(scale_mode, "global")) {
        scale_margin <- integer()
    } else if (identical(scale_mode, "product")) {
        scale_margin <- which(can_scale)
        if (length(scale_margin) == 0L) {
            stop("separable() scale = product() needs at least one ",
                 "scale-capable margin.")
        }
    } else if (identical(scale_mode, "selected_product")) {
        scale_spec <- scale$margins
        scale_margin <- match(.sep_margin_key(scale_spec),
                              .sep_margin_key(margins))
        if (anyNA(scale_margin)) {
            i <- which(is.na(scale_margin))[[1]]
            stop("separable() scale margin ", scale_spec$struc[i], "(",
                 scale_spec$var[i], ") must match one of the specified ",
                 "margins.")
        }
        scale_ok <- vapply(regs[scale_margin],
                            function(x) x$scale$can_scale, logical(1))
        if (!all(scale_ok)) {
            i <- which(!scale_ok)[[1]]
            stop("separable() scale = ", scale_spec$struc[i], "(",
                 scale_spec$var[i], ") selects a correlation-only margin. ",
                 "Use a scale-capable margin or scale = global().")
        }
        if (anyDuplicated(scale_margin)) {
            stop("separable() scale margins must be unique.")
        }
    } else if (identical(scale_mode, "none")) {
        scale_margin <- integer()
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
    supported <- names(.sep_margin_registry)[
        vapply(.sep_margin_registry, function(x) {
            x$builder %in% names(.sep_builder_kind_code)
        }, logical(1))
    ]
    stop("separable() frontend parsed ", .sep_margin_label(margins),
         ", but the backend currently only evaluates products among ",
         paste0(supported, "()", collapse = ", "), ".")
}

.sep_distance_info <- function(regs, spec, dims) {
    needs_dist <- vapply(regs, function(x) x$metadata$needs_dist, logical(1))
    starts <- rep.int(-1L, length(dims))
    dists <- numeric()
    if (!any(needs_dist)) return(list(starts = starts, dists = dists))

    if (is.null(spec$margin_cnms)) {
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
        starts[[i]] <- length(dists)
        dists <- c(dists, as.vector(as.matrix(stats::dist(coords))))
    }
    list(starts = starts, dists = dists)
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

    if (length(rows)) do.call(rbind, rows) else
        data.frame(margin = integer(), kind = integer(),
                   start = integer(), length = integer())
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

    theta_layout <- .sep_theta_block_layout(regs, dims, scale_info)
    ntheta <- sum(theta_layout$length)
    builder_kind <- vapply(regs, `[[`, character(1), "builder")
    scale_kind <- vapply(regs, function(x) x$scale$kind, character(1))
    distance_info <- .sep_distance_info(regs, spec, dims)
    cov_info <- .sep_fixed_cov_info(regs, spec, dims)

    list(
        dims = dims,
        codes = as.integer(vapply(strucs, function(z) .valid_covstruct[[z]], numeric(1))),
        builder_kinds = as.integer(.sep_builder_kind_code[builder_kind]),
        scale_kinds = as.integer(.sep_scale_kind_code[scale_kind]),
        dispatch = as.integer(.sep_dispatch_code[dispatch]),
        scale_mode = scale_info$mode_code,
        scale_spec = scale_info$spec,
        theta_block_margins = theta_layout$margin,
        theta_block_kinds = theta_layout$kind,
        theta_block_starts = theta_layout$start,
        theta_block_lengths = theta_layout$length,
        dist_starts = distance_info$starts,
        dists = distance_info$dists,
        fixed_cov_starts = cov_info$starts,
        fixed_covs = cov_info$covs,
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
    if (!all(margins$struc %in% names(.sep_margin_registry))) {
        bad <- unique(margins$struc[!margins$struc %in% names(.sep_margin_registry)])
        stop("Unsupported separable() margin: ", paste(bad, collapse = ", "))
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

.sep_check_cov_payload <- function(x, cnms, expr, struc, env) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("separable() ", struc, "() margin expects a numeric matrix ",
             "extra argument.", call. = FALSE)
    }
    if (nrow(x) != ncol(x)) {
        stop("separable() ", struc, "() margin covariance matrix must be square.",
             call. = FALSE)
    }
    if (nrow(x) != length(cnms)) {
        stop("separable() ", struc, "() margin covariance matrix has ",
             "dimension ", nrow(x), ", but the margin has ", length(cnms),
             " columns.", call. = FALSE)
    }
    if (identical(struc, "propto")) {
        rn <- rownames(x)
        cn <- colnames(x)
        if (is.null(rn) && is.null(cn)) {
            stop("row or column names of propto matrix are required",
                 call. = FALSE)
        }
        if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
            stop("row and column names of propto matrix do not match",
                 call. = FALSE)
        }
        mat_names <- if (is.null(cn)) rn else cn
        if (!identical(mat_names, cnms)) {
            labs <- attr(stats::terms(.sep_margin_formula(expr, env)),
                         "term.labels")
            prefixed <- paste0(labs, mat_names)
            if (!identical(prefixed, cnms)) {
                stop("column or row names of the propto matrix do not match ",
                     "the separable() margin columns.", call. = FALSE)
            }
        }
    }
    x
}

.sep_check_payload <- function(payload, value, cnms, expr, struc, env) {
    switch(payload$kind,
           cov_matrix = .sep_check_cov_payload(value, cnms, expr, struc, env),
           stop("Unsupported separable() payload kind: ", payload$kind,
                call. = FALSE))
}

.sep_resolve_margin_payloads <- function(margins, Xlist, fr, env) {
    payloads <- vector("list", nrow(margins))
    for (i in seq_len(nrow(margins))) {
        reg <- .sep_margin_registry[[margins$struc[i]]]
        specs <- reg$extra$payloads
        if (!length(specs)) next
        payloads[[i]] <- setNames(lapply(specs, function(payload) {
            value <- .sep_eval_extra(margins$extra[[i]][[payload$arg]], fr, env)
            .sep_check_payload(payload, value, colnames(Xlist[[i]]),
                               margins$expr[[i]], margins$struc[i], env)
        }), vapply(specs, `[[`, character(1), "kind"))
    }
    payloads
}

.sep_fixed_cov_info <- function(regs, spec, dims) {
    starts <- rep.int(-1L, length(dims))
    covs <- numeric()
    for (i in seq_along(regs)) {
        if (!identical(regs[[i]]$builder, "fixed_cov")) next
        cov <- spec$margin_payloads[[i]][["cov_matrix"]]
        if (is.null(cov) || nrow(cov) != dims[[i]] || ncol(cov) != dims[[i]]) {
            stop("separable() fixed covariance margin metadata does not match ",
                 "the product design.")
        }
        starts[[i]] <- length(covs)
        covs <- c(covs, as.vector(cov))
    }
    list(starts = starts, covs = covs)
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
