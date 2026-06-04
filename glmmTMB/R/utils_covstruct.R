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
##' ## member/time helper for product covariance structures
##' membertime(factor(c("A", "B", "A", "B")), c(1, 1, 2, 2))
##' mt(factor(c("A", "B", "A", "B")), c(1, 1, 2, 2))
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
##' @param member member or role coordinate for member-by-time covariance structures
##' @param time discrete time coordinate for member-by-time covariance structures
##' @details \code{membertime} is a convenience wrapper for creating two-dimensional
##' coordinate factors for product covariance structures such as
##' \code{homcsxar1} and \code{unxar1}. It stores the first coordinate as the
##' member/role index and the second coordinate as the time index. The first
##' coordinate may have two or more levels; dyads are the common motivating
##' case, but larger groups are allowed. If
##' \code{member} is a factor, its levels are retained for compact
##' \code{print} and \code{summary} output. For \code{homcsxar1}, member
##' labels are exchangeable and may be arbitrary, but the first coordinate
##' must identify the same individual within each group across time. If the
##' mean model uses an \code{Idiff} coding for indistinguishable dyads, use
##' the same stable member assignment in \code{membertime}.
##' @export
membertime <- function(member, time) {
    coord <- function(x, time = FALSE) {
        if (is.factor(x)) {
            if (time) {
                xx <- suppressWarnings(as.numeric(as.character(x)))
                if (!anyNA(xx)) return(xx)
            }
            return(as.numeric(x))
        }
        x
    }
    member_coord <- coord(member)
    ans <- numFactor(member_coord, coord(time, time = TRUE))
    member_values <- sort(unique(member_coord))
    member_labels <- if (is.factor(member)) {
        levels(member)[member_values]
    } else {
        as.character(member_values)
    }
    names(member_labels) <- as.character(member_values)
    attr(ans, "memberLabels") <- member_labels
    ans
}

##' @rdname numFactor
##' @export
mt <- membertime

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
