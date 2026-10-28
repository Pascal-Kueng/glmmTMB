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

## How a kron() term is fitted. Each margin starts as a spec (a list), and
## four functions handle it in turn:
## 1. rewriteKron(), below, makes the spec from the formula: struc (e.g.
##    "us"), lhs (e.g. 1 + x) and arg (e.g. propto's matrix, or NA).
## 2. getXReTrms() adds cnms (the margin's column names) and builds the
##    term's Z as the Khatri-Rao product of the margins' model matrices.
## 3. getReStruc() turns each spec into a covariance term by calling itself,
##    and adds theta (start values) and free (which of them are estimated).
##    Later margins get their scale fixed here.
## 4. kron_nll() fills the estimated values into each margin's theta, then
##    combines the margins' densities with RTMB::dseparable() and applies
##    any rr loadings.

## Rewrite each kron(A(lhs1) %x% B(lhs2, arg) | g) in a formula into the
## internal form
##   kron(lhs1 + lhs2 | g, list(list(struc = "A", lhs = quote(lhs1), arg = NA),
##                              list(struc = "B", lhs = quote(lhs2), arg = arg)))
## parsed like rr(0 + f | g, 2). The first argument only puts the margins'
## variables into the model frame: getXReTrms() builds Z from the margins.
rewriteKron <- function(x) {
    ## Rewrite only the one-argument user syntax, not the internal form (as
    ## returned by formula(fit, component = "cond")). Otherwise, look for
    ## kron() terms in the arguments.
    if (!identical(x[[1]], quote(kron)) || length(x) != 2) {
        for (i in seq_along(x)[-1]) if (is.call(x[[i]])) x[[i]] <- rewriteKron(x[[i]])
        return(x)
    }
    bar <- x[[2]] # A(lhs1) %x% B(lhs2, arg) | g
    ## Require "| g" (a "|| g" is read as "| g", as for other structures)
    if (!is.call(bar) || !deparse1(bar[[1]]) %in% c("|", "||"))
        stop("kron() needs a grouping factor, as in ",
             "kron(us(0 + f1) %x% ar1(0 + f2) | g)", call. = FALSE)
    ## A(lhs1) %x% B(lhs2, arg) %x% ...  ->  list(A(lhs1), B(lhs2, arg), ...)
    flatten <- function(e) {
        if (!is.call(e) || !identical(e[[1]], quote(`%x%`))) return(list(e))
        c(flatten(e[[2]]), flatten(e[[3]]))
    }
    margins <- flatten(bar[[2]])
    if (!all(vapply(margins, function(e) is.call(e) && length(e) %in% 2:3, NA)))
        stop("kron() margins must look like us(1 + x) or propto(0 + f, C)",
             call. = FALSE)
    lhs <- Reduce(function(a, b) call("+", a, b), lapply(margins, `[[`, 2))
    if (any(c("|", "||") %in% all.names(lhs)))
        stop("kron() margins take no grouping bar: use one outer bar, as in ",
             "kron(us(0 + f1) %x% ar1(0 + f2) | g)", call. = FALSE)
    ## A(lhs, arg)  ->  list(struc = "A", lhs = quote(lhs), arg = arg)
    marginSpec <- function(e)
        call("list", struc = deparse(e[[1]]), lhs = call("quote", e[[2]]),
             arg = if (length(e) == 3) e[[3]] else NA)
    spec <- as.call(c(quote(list), lapply(margins, marginSpec)))
    ## Split nested groups here (g/h -> h:g and g, as for us(0 + f | g/h)):
    ## splitForm() does not split them in terms with a second argument
    groups <- reformulas::expandGrpVar(bar[[3]])
    ## "| 1" expands to no groups, which would drop the term silently
    if (!length(groups))
        stop("kron() needs a grouping factor: for a single group, use a ",
             "factor with one level instead of | 1", call. = FALSE)
    reformulas::sumTerms(lapply(groups, function(g)
        call("kron", call("|", lhs, g), spec)))
}

