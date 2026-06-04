## DHARMa diagnostics for latent separable member-by-time AR(1) terms.

get_xar1_innovations <- function(fit, component = "cond", term = NULL) {
    if (!inherits(fit, "glmmTMB")) {
        stop("'fit' must be a glmmTMB model.")
    }
    component <- match.arg(component, c("cond", "zi", "disp"))
    vc_list <- VarCorr(fit)[[component]]
    is_xar1 <- vapply(
        vc_list,
        function(x) any(class(x) %in% c("vcmat_homcsxar1", "vcmat_unxar1")),
        logical(1)
    )
    if (!any(is_xar1)) {
        stop("No homcsxar1 or unxar1 term found in the fitted model.")
    }
    xar1_names <- names(vc_list)[is_xar1]
    vc_name <- if (is.null(term)) {
        if (length(xar1_names) > 1) {
            warning("Multiple homcsxar1/unxar1 terms found; using the first. ",
                    "Specify 'term' to choose a term explicitly.")
        }
        xar1_names[1]
    } else if (is.numeric(term)) {
        xar1_names[term]
    } else if (term %in% names(vc_list)) {
        term
    } else {
        stop("'term' must be NULL, a numeric index, or one of: ",
             paste(xar1_names, collapse = ", "))
    }
    if (is.na(vc_name) || !vc_name %in% xar1_names) {
        stop("'term' does not select a homcsxar1 or unxar1 term.")
    }

    vc <- vc_list[[vc_name]]
    sd <- attr(vc, "stddev")
    cor <- attr(vc, "correlation")
    coords <- parseNumLevels(names(sd))
    if (ncol(coords) != 2) {
        stop("xar1 terms require two-dimensional member-time coordinates.")
    }
    members <- sort(unique(coords[, 1]))
    times <- sort(unique(coords[, 2]))
    if (length(times) < 2) {
        stop("At least two time points are required to compute AR(1) innovations.")
    }

    pos <- function(member, time) {
        which(coords[, 1] == member & coords[, 2] == time)[1]
    }
    member_pos <- vapply(members, pos, integer(1), time = times[1])
    phi <- cor[pos(members[1], times[1]), pos(members[1], times[2])]
    Sigma_member <- outer(sd[member_pos], sd[member_pos]) *
        cor[member_pos, member_pos, drop = FALSE]
    Sigma_innov <- (1 - phi^2) * Sigma_member

    re_group <- sub("\\.[[:digit:]]+$", "", vc_name)
    re <- ranef(fit, condVar = FALSE)[[component]][[re_group]]
    cols <- names(sd)
    if (is.null(re) || !all(cols %in% colnames(re))) {
        stop("Could not match the covariance term columns to ranef() output.")
    }

    u <- as.matrix(re[, cols, drop = FALSE])
    u_array <- array(
        NA_real_,
        dim = c(nrow(u), length(times), length(members)),
        dimnames = list(rownames(u), as.character(times), as.character(members))
    )
    for (j in seq_along(cols)) {
        u_array[, match(coords[j, 2], times), match(coords[j, 1], members)] <- u[, j]
    }

    innov <- do.call(
        rbind,
        lapply(seq_len(dim(u_array)[1]), function(g) {
            cur <- matrix(u_array[g, -1, ], ncol = length(members))
            prev <- matrix(u_array[g, -length(times), ], ncol = length(members))
            cur - phi * prev
        })
    )
    z <- t(forwardsolve(t(chol(Sigma_innov)), t(innov)))
    colnames(innov) <- colnames(z) <- colnames(Sigma_member) <-
        rownames(Sigma_member) <- as.character(members)

    list(
        term = vc_name,
        phi = phi,
        Sigma_member = Sigma_member,
        Sigma_innov = Sigma_innov,
        innovations = innov,
        standardized = z
    )
}

##' Wrap separable AR(1) latent innovations for DHARMa diagnostics
##'
##' \code{dharma_xar1} returns a lightweight wrapper around a fitted
##' \code{glmmTMB} model with a \code{homcsxar1} or \code{unxar1} random-effect
##' term. The wrapper lets \code{DHARMa::simulateResiduals()} diagnose the
##' standardized one-step latent AR(1) innovations rather than the observed
##' response. This is useful for Gaussian models fitted with
##' \code{dispformula = ~0}, where the lowest-level residual process is modeled
##' as a structured latent random effect.
##'
##' @param fit a fitted \code{glmmTMB} model containing a \code{homcsxar1} or
##' \code{unxar1} term
##' @param component model component; currently defaults to the conditional
##' component
##' @param term optional \code{homcsxar1}/\code{unxar1} term name or numeric
##' index if the model contains more than one such term
##' @return an object suitable for \code{DHARMa::simulateResiduals()}
##' @details The resulting diagnostics are based on fitted conditional modes of
##' the latent random effects. They are useful model checks, but should not be
##' interpreted as ordinary response-scale residual diagnostics.
##' @examples
##' \dontrun{
##' res <- DHARMa::simulateResiduals(dharma_xar1(fit))
##' plot(res)
##' }
##' @export
dharma_xar1 <- function(fit, component = "cond", term = NULL) {
    x <- get_xar1_innovations(fit, component = component, term = term)
    attr(fit, "dharma_xar1") <- x
    fit
}

xar1_dharma_data <- function(object) {
    attr(object, "dharma_xar1")
}

xar1_dharma_response <- function(object) {
    as.vector(xar1_dharma_data(object)$standardized)
}

xar1_dharma_simulate <- function(object, nsim = 1, seed = NULL) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        runif(1)
    }
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ret <- as.data.frame(replicate(nsim, rnorm(length(xar1_dharma_response(object)))))
    names(ret) <- paste0("sim_", seq_len(nsim))
    attr(ret, "seed") <- RNGstate
    ret
}
