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

    struc <- character(length(margins))
    expr <- vector("list", length(margins))
    margin_names <- character(length(margins))

    for (i in seq_along(margins)) {
        margin <- margins[[i]]
        arg_names <- names(margin)[-1L]
        has_named_arg <- length(arg_names) > 0L && any(nzchar(arg_names))
        if (!is.call(margin) || length(margin) != 2L || has_named_arg) {
            stop("kron() margins must look like us(1 + x) or ",
                 "ar1(0 + time)", call. = FALSE)
        }
        struc[[i]] <- .kron_call_name(margin)
        if (!struc[[i]] %in% .kron_supported) {
            stop("unsupported kron() margin: ", struc[[i]], call. = FALSE)
        }
        expr[[i]] <- margin[[2L]]
        margin_names[[i]] <- deparse1(expr[[i]], collapse = "")
    }

    list(struc = struc, expr = expr, margin_names = margin_names)
}
