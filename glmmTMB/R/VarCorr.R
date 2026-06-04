## returns a true family() object iff one was given
## to glmmTMB() in the first place ....
##' @importFrom stats family
##' @export
family.glmmTMB <- function(object, ...) {
    if (!is.null(xar1_dharma_data(object))) return(gaussian())
    object$modelInfo$family
}

## don't quite know why this (rather than just ...$parList()) is
## necessary -- used in ranef.glmmTMB and sigma.glmmTMB
getParList <- function(object) {
    object$obj$env$parList(object$fit$par, object$fit$parfull)
}


##' Extract residual standard deviation or dispersion parameter
##'
##' For Gaussian models, \code{sigma} returns the value of the residual
##' standard deviation; for other families, it returns the
##' dispersion parameter, \emph{however it is defined for that
##' particular family}. See details for each family below.
##'
##' @details
##'  The value returned varies by family:
##'  \describe{
##'      \item{gaussian}{returns the \emph{maximum likelihood} estimate
##'          of the standard deviation (i.e., smaller than the results of
##'                                  \code{sigma(lm(...))} by a factor of (n-1)/n)}
##'      \item{nbinom1}{returns a dispersion parameter
##'          (usually denoted \eqn{\alpha}{alpha} as in Hardin and Hilbe (2007)):
##'          such that the variance equals \eqn{\mu(1+\alpha)}{mu(1+alpha)}.}
##'      \item{nbinom2}{returns a dispersion parameter
##'          (usually denoted \eqn{\theta}{theta} or \eqn{k}); in contrast to
##'          most other families, larger \eqn{\theta}{theta} corresponds to a \emph{lower}
##'          variance which is \eqn{\mu(1+\mu/\theta)}{mu(1+mu/theta)}.}
##'      \item{Gamma}{Internally, glmmTMB fits Gamma responses by fitting a mean
##'          and a shape parameter; sigma is estimated as (1/sqrt(shape)),
##'          which will typically be close (but not identical to) that estimated
##'          by \code{stats:::sigma.default}, which uses sqrt(deviance/df.residual)}
##'      \item{beta}{returns the value of \eqn{\phi}{phi},
##'          where the conditional variance is \eqn{\mu(1-\mu)/(1+\phi)}{mu*(1-mu)/(1+phi)}
##'          (i.e., increasing \eqn{\phi}{phi} decreases the variance.)
##'          This parameterization follows Ferrari and Cribari-Neto (2004)
##'          (and the \code{betareg} package):}
##'      \item{betabinomial}{This family uses the same parameterization (governing
##'           the Beta distribution that underlies the binomial probabilities) as \code{beta}.}
##'      \item{genpois}{returns the index of dispersion \eqn{\phi^2}{phi^2},
##'           where the variance is \eqn{\mu\phi^2}{mu*phi^2} (Consul & Famoye 1992)}
##'      \item{compois}{returns the value of \eqn{1/\nu}{1/nu};
##'           when \eqn{\nu=1}{nu=1}, compois is equivalent to the Poisson distribution.
##'           There is no closed form equation for the variance, but
##'           it is approximately underdispersed when \eqn{1/\nu <1}{1/nu <1}
##'           and approximately overdispersed when \eqn{1/\nu >1}{1/nu>1}.
##'           In this implementation, \eqn{\mu}{mu} is exactly equal to the mean (Huang 2017), which
##'           differs from the COMPoissonReg package (Sellers & Lotze 2015).}
##'      \item{tweedie}{returns the value of \eqn{\phi}{phi},
##'           where the variance is \eqn{\phi\mu^p}{phi*mu^p}.
##'           The value of \eqn{p} can be extracted using \code{family_params}
##'      }
##'      \item{ordbeta}{see details for \code{beta}}
##' }
##'
##'  The most commonly used GLM families
##'  (\code{binomial}, \code{poisson}) have fixed dispersion parameters which are
##'  internally ignored.
##'
##' @references
##' \itemize{
##' \item Consul PC, and Famoye F (1992). "Generalized Poisson regression model. Communications in Statistics: Theory and Methods" 21:89–109.
##' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
##' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions." Stata press.
##' \item Huang A (2017). "Mean-parametrized Conway–Maxwell–Poisson regression models for dispersed counts. " \emph{Statistical Modelling} 17(6), 1-22.
##' \item Sellers K & Lotze T (2015). "COMPoissonReg: Conway-Maxwell Poisson (COM-Poisson) Regression". R package version 0.3.5. https://CRAN.R-project.org/package=COMPoissonReg
##' }
##' @aliases sigma
##' @param object a \dQuote{glmmTMB} fitted object
##' @param \dots (ignored; for method compatibility)
## Import generic and re-export
## note the following line is hacked in Makefile/namespace-update to ...
## if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)
## also see <https://github.com/klutometis/roxygen/issues/371>
##' @rawNamespace if(getRversion()>='3.3.0') { importFrom(stats, sigma) } else { importFrom(lme4,sigma) }
##  n.b. REQUIRES roxygen2 >= 5.0
## @importFrom lme4 sigma
##' @export sigma
##' @method sigma glmmTMB
##' @export
sigma.glmmTMB <- function(object, ...) {
    pl <- getParList(object)
    ff <- object$modelInfo$family$family
    if (!usesDispersion(ff)) return(1.)
    if (length(pl$betadisp)>1) return(NA)
    switch(family(object)$family,
           Gamma=exp(-0.5*pl$betadisp),
           exp(pl$betadisp))
}

##' @noRd
##' @param cor list of correlation matrices/structures
##' @param sd list of sd vectors
##' @param cnms vector of term names
##' @param sc scale parameter (sigma)
##' @param bc vector of block codes
##' @param useSc use scale parameter?
mkVC <- function(cor, sd, cnms, sc, bc, useSc) {
    stopifnot(length(cnms) == (nc <- length(cor)),  nc == length(sd),
              is.list(cnms), is.list(cor), is.list(sd),
              is.character(nnms <- names(cnms)), nzchar(nnms))

    ## FIXME: do we want this?  Maybe not.
    ## Potential problem: the names of the elements of the VarCorr() list
    ##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
    ##  Subject terms, so the names are "Subject", "Subject".  The print method
    ##  for VarCorrs handles this just fine, but it's a little awkward if we
    ##  want to dig out elements of the VarCorr list ... ???
    if (anyDuplicated(nnms))
        nnms <- make.names(nnms, unique = TRUE)
    ##
    ## cov :=  F(sd, cor) :
    do1cov <- function(sd, cor, n = length(sd)) {
        sd * cor * rep(sd, each = n)
    }
    docov <- function(sd, cor, nm, bc) {
        maxdim <- max(length(sd), nrow(cor))
        ## extend sd if necessary
        if (length(sd)==1 && maxdim > 1) {
            sd <- rep(sd, maxdim)
        }
        ## truncate names if necessary
        if (length(nm) > maxdim) {
            nm <- nm[seq(maxdim)]
        }
        ## what do we do with models without fully filled in cor structures?
        ## (diag, ar1, cs, ...)
        diagmodel <- identical(dim(cor),c(0L,0L))
        if (diagmodel) cor <- diag(length(sd))
        cov <- matrix(numeric(0))
        if (identical(c(cor), NaN)) {
            cov <- NA
        } else {
            if (length(sd) == 1 || length(sd) == nrow(cor)) {
                cov <- do1cov(sd, cor)
            }
        }
        ## may not want full names (e.g. ou/ar1 models)
        names(sd) <- nm[seq_along(sd)]
        ## skip naming if null cov/cor matrix
        if (!identical(cov, NA) && !is.null(cov) && nrow(cov) > 0) {
            dimnames(cov) <- dimnames(cor) <- list(nm,nm)
        }
        structure(cov, stddev=sd, correlation=cor)
    }
    ss <- setNames(mapply(docov, sd, cor, cnms, bc, SIMPLIFY=FALSE),nnms)
    ## ONLY first element -- otherwise breaks formatVC
    ## FIXME: do we want a message/warning here, or elsewhere,
    ##   when the 'Residual' var parameters are truncated?
    attr(ss,"sc") <- sc[1]
    attr(ss,"useSc") <- useSc
    ss
}

xar1_block_names <- c("unxar1", "homcsxar1")

is_xar1_vcmat <- function(x) {
    bc <- attr(x, "blockCode")
    !is.null(bc) && names(bc)[1] %in% xar1_block_names
}

xar1_term_label <- function(x) {
    term <- sub("\\([^()]*\\)$", "", x)
    helper <- regexec("^(?:(?:glmmTMB::)?)(?:membertime|mt)\\((.*)\\)$", term)
    parts <- regmatches(term, helper)[[1]]
    if (length(parts) == 2) {
        args <- trimws(strsplit(parts[2], ",", fixed = TRUE)[[1]])
        args <- trimws(sub("^[[:alnum:]_.]+[[:space:]]*=[[:space:]]*", "", args))
        if (length(args) >= 2) return(paste(args[1], "x", args[2]))
    }
    term
}

##' Extract variance and correlation components
##'
##' @aliases VarCorr
##' @param x a fitted \code{glmmTMB} model
##' @param sigma residual standard deviation (usually set automatically from internal information)
##' @param ... extra arguments (for consistency with generic method)
##' @importFrom nlme VarCorr
## and re-export the generic:
##' @export VarCorr
##' @export
##' @examples
##' ## Comparing variance-covariance matrix with manual computation
##' data("sleepstudy",package="lme4")
##' fm4 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
##' VarCorr(fm4)[[c("cond","Subject")]]
##' ## hand calculation
##' pars <- getME(fm4,"theta")
##' ## construct cholesky factor
##' L <- diag(2)
##' L[lower.tri(L)] <- pars[-(1:2)]
##' C <- crossprod(L)
##' diag(C) <- 1
##' sdvec <- exp(pars[1:2])
##' (V <- outer(sdvec,sdvec) * C)
##' @details For an unstructured variance-covariance matrix, the internal parameters
##' are structured as follows: the first n parameters are the log-standard-deviations,
##' while the remaining n(n-1)/2 parameters are the elements of the Cholesky factor
##' of the correlation matrix, filled in column-wise order
##' (see the \href{http://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html}{TMB documentation}
##' for further details).
##' @keywords internal
VarCorr.glmmTMB <- function(x, sigma = 1, ... )
{
    ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    ## FIXME:: do we need 'sigma' any more (now that nlme generic
    ##         doesn't have it?)
    check_dots(..., .action = "warning")
    stopifnot(is.numeric(sigma), length(sigma) == 1)
    xrep <- x$obj$env$report(x$fit$parfull)
    reT <- x$modelInfo$reTrms
    reS <- x$modelInfo$reStruc
    familyStr <- family(x)$family
    useSc <- if (missing(sigma)) {
                 ## *only* report residual variance for Gaussian family ...
                 ## *not* usesDispersion(familyStr)
                 sigma <- sigma(x)
                 familyStr=="gaussian" && !zeroDisp(x)
             } else TRUE
    comp_nms <- c("", "zi", "disp")
    comp_nms2 <- c("cond", "zi", "disp")
    corr_list <- vector(mode="list", length=length(comp_nms2))
    names(corr_list) <- comp_nms2
    get_member_labels <- function(term) {
        nm <- sub("\\([^()]*\\)$", "", term[1])
        if (is.null(x$frame)) return(NULL)
        if (!nm %in% names(x$frame)) return(NULL)
        attr(x$frame[[nm]], "memberLabels")
    }

    for (i in seq_along(comp_nms)) {
        restruc <- reS[[paste0(comp_nms2[i],  "ReStruc")]]
        ## lapply() rather than [vs]apply, don't want to lose names
        bcvec <- lapply(restruc, function(x) x[["blockCode"]])
        if(length(cn <- reT[[comp_nms2[i]]]$cnms)) {
            vc <- mkVC(cor = xrep[[paste0("corr", comp_nms[i])]],
                       sd  = xrep[[paste0("sd", comp_nms[i])]],
                       cnms = cn,
                       sc = sigma,
                       bc = bcvec,
                       useSc = useSc)
            for (j in seq_along(vc)) {
                attr(vc[[j]],"blockCode") <- bcvec[[j]]
                if (names(bcvec[[j]]) %in% xar1_block_names) {
                    attr(vc[[j]], "memberLabels") <- get_member_labels(cn[[j]])
                }
                class(vc[[j]]) <- c(paste0("vcmat_", names(bcvec[[j]])), class(vc[[j]]))
            }
            corr_list[[comp_nms2[[i]]]] <- vc
        }
    }
    structure(corr_list,
	      sc = usesDispersion(familyStr), ## 'useScale'
	      class = c("VarCorr.glmmTMB"))
}

##' Printing The Variance and Correlation Parameters of a \code{glmmTMB}
##' @method print VarCorr.glmmTMB
##' @importFrom reformulas formatVC
##' @export
##  document as it is a method with "surprising arguments":
##' @param x a result of \code{\link{VarCorr}(<glmmTMB>)}.
##' @param digits number of significant digits to use.
##' @param comp a string specifying the component to format and print.
##' @param formatter a \code{\link{function}}.
##' @param maxdim maximum number of SDs/correlations to print
##' @param ... optional further arguments, passed to \code{\link{print.default}} (for printing a character matrix)
print.VarCorr.glmmTMB <- function(x, digits = max(3, getOption("digits") - 2),
				  comp = "Std.Dev.", formatter = format, maxdim = 10, ...)
{
  for (cc in names(x)) {
    if(!is.null(x[[cc]])) {
      cat(sprintf("\n%s:\n", cNames[[cc]]))
      print(formatVC(x[[cc]],
                     digits = digits, comp = comp, formatter = formatter, maxdim = maxdim),
            quote = FALSE, ...)
    }
  }
  invisible(x)
}

formatVC.glmmTMB <- function(x, digits = max(3, getOption("digits") - 2),
                             comp = "Std.Dev.", formatter = format,
                             maxdim = 10, compactXAr1 = FALSE) {
    if (!compactXAr1 || !any(vapply(x, is_xar1_vcmat, logical(1)))) {
        return(formatVC(x, digits = digits, comp = comp,
                        formatter = formatter, maxdim = maxdim))
    }
    ff <- function(z) formatter(z, digits = digits)
    empty_row <- function() setNames(rep("", length(c("Groups", "Name", comp, "Corr"))),
                                     c("Groups", "Name", comp, "Corr"))
    fmt_comp <- function(sd) {
        ans <- setNames(rep("", length(comp)), comp)
        if ("Variance" %in% comp) ans["Variance"] <- ff(sd^2)
        if ("Std.Dev." %in% comp) ans["Std.Dev."] <- ff(sd)
        ans
    }
    pos <- function(coords, member, time) {
        which(coords[, 1] == member & coords[, 2] == time)[1]
    }
    format_xar1 <- function(z, group) {
        bc <- names(attr(z, "blockCode"))[1]
        sd <- attr(z, "stddev")
        cor <- attr(z, "correlation")
        member_labels <- attr(z, "memberLabels")
        coords <- parseNumLevels(names(sd))
        members <- sort(unique(coords[, 1]))
        times <- sort(unique(coords[, 2]))
        term <- xar1_term_label(names(sd)[1])
        first_time <- times[1]
        first_member <- members[1]
        row_list <- list()
        member_name <- function(m) {
            if (!is.null(member_labels)) {
                lab <- member_labels[as.character(m)]
                if (!is.na(lab)) return(unname(lab))
            }
            as.character(m)
        }

        member_pos <- vapply(members, function(m) pos(coords, m, first_time), integer(1))
        if (bc == "homcsxar1") {
            rr <- empty_row()
            rr["Groups"] <- group
            rr["Name"] <- paste0(term, " member")
            rr[comp] <- fmt_comp(sd[member_pos[1]])
            if (length(members) > 1) {
                rho <- cor[member_pos[1], member_pos[2]]
                rr["Corr"] <- paste(ff(rho), "(homcs)")
            }
            row_list[[length(row_list) + 1]] <- rr
        } else {
            for (i in seq_along(members)) {
                rr <- empty_row()
                if (i == 1) rr["Groups"] <- group
                rr["Name"] <- paste0(term, " ", member_name(members[i]))
                rr[comp] <- fmt_comp(sd[member_pos[i]])
                if (i > 1) {
                    rr["Corr"] <- paste(ff(cor[member_pos[i], member_pos[seq_len(i - 1)]]),
                                        collapse = " ")
                }
                row_list[[length(row_list) + 1]] <- rr
            }
        }

        if (length(times) > 1) {
            phi <- cor[pos(coords, first_member, times[1]),
                       pos(coords, first_member, times[2])]
            rr <- empty_row()
            rr["Name"] <- paste0(term, " time lag 1")
            rr["Corr"] <- paste(ff(phi), "(ar1)")
            row_list[[length(row_list) + 1]] <- rr

            if (length(members) > 1) {
                pairs <- utils::combn(members, 2)
                for (j in seq_len(ncol(pairs))) {
                    rr <- empty_row()
                    rr["Name"] <- paste0(term, " members ",
                                         member_name(pairs[1, j]), ":",
                                         member_name(pairs[2, j]),
                                         " lag 1")
                    rr["Corr"] <- paste(
                        ff(cor[pos(coords, pairs[1, j], times[1]),
                               pos(coords, pairs[2, j], times[2])]),
                        "(xar1)")
                    row_list[[length(row_list) + 1]] <- rr
                }
            }
        }
        do.call(rbind, row_list)
    }
    rows <- list()
    for (i in seq_along(x)) {
        xi <- x[[i]]
        if (is_xar1_vcmat(xi)) {
            rows[[length(rows) + 1]] <- format_xar1(xi, names(x)[i])
        } else {
            xx <- x[i]
            attr(xx, "useSc") <- FALSE
            rows[[length(rows) + 1]] <- formatVC(xx, digits = digits,
                                                 comp = comp,
                                                 formatter = formatter,
                                                 maxdim = maxdim)
        }
    }
    if (isTRUE(attr(x, "useSc"))) {
        rr <- empty_row()
        rr["Groups"] <- "Residual"
        rr[comp] <- fmt_comp(attr(x, "sc"))
        rows[[length(rows) + 1]] <- matrix(rr, nrow = 1,
                                           dimnames = list(NULL, names(rr)))
    }
    max_cols <- max(vapply(rows, ncol, integer(1)))
    col_names <- colnames(rows[[which.max(vapply(rows, ncol, integer(1)))]])
    pad_cols <- function(m) {
        if (ncol(m) < max_cols) {
            extra <- matrix("", nrow = nrow(m), ncol = max_cols - ncol(m))
            m <- cbind(m, extra)
        }
        colnames(m) <- col_names
        m
    }
    rows <- lapply(rows, pad_cols)
    ans <- do.call(rbind, rows)
    rownames(ans) <- rep("", nrow(ans))
    ans
}
