## DHARMa diagnostics for Gaussian models with separable member-by-time AR(1) terms.

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

get_xar1_term_name <- function(fit, term = NULL) {
    vc_list <- VarCorr(fit)$cond
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
    vc_name
}

chol_with_jitter <- function(x) {
    R <- try(chol(x), silent = TRUE)
    if (!inherits(R, "try-error")) return(R)
    eps <- sqrt(.Machine$double.eps) * max(1, max(diag(x)))
    chol(x + diag(eps, nrow(x)))
}

get_xar1_marginal_residuals <- function(fit, term = NULL) {
    if (!inherits(fit, "glmmTMB")) {
        stop("'fit' must be a glmmTMB model.")
    }
    fam <- family(fit)
    if (!(fam$family == "gaussian" && fam$link == "identity")) {
        stop("dharma_xar1() currently provides marginal normalized residuals ",
             "only for Gaussian identity-link models. For count or other ",
             "non-Gaussian models, use DHARMa::simulateResiduals(fit, ",
             "simulateREs = \"conditional\") for response-scale diagnostics.")
    }
    vc_name <- get_xar1_term_name(fit, term = term)

    reTrms <- fit$modelInfo$reTrms$cond
    if (is.null(reTrms$cnms)) {
        stop("No conditional random effects found.")
    }
    flist <- reTrms$flist
    if (length(flist) != 1 || !all(attr(flist, "assign") == 1)) {
        stop("dharma_xar1() currently supports Gaussian models whose ",
             "conditional random-effect terms all use a single grouping factor.")
    }

    group <- flist[[1]]
    nlev <- nlevels(group)
    Z <- getME(fit, "Z")
    reStruc <- fit$modelInfo$reStruc$condReStruc
    Gp <- cumsum(c(0, vapply(reStruc, function(x) x$blockReps * x$blockSize,
                             numeric(1))))
    if (length(Gp) != length(reStruc) + 1) {
        stop("Could not align random-effect structure with the Z matrix.")
    }
    if (!all(vapply(reStruc, `[[`, numeric(1), "blockReps") == nlev)) {
        stop("dharma_xar1() currently requires one random-effect block per ",
             "level of the grouping factor for every conditional term.")
    }

    y <- model.response(fit$frame)
    if (!is.numeric(y) || !is.null(dim(y))) {
        stop("dharma_xar1() requires a numeric Gaussian response.")
    }
    X <- getME(fit, "X")
    beta <- getParList(fit)$beta
    mu <- as.vector(X %*% beta)
    e <- y - mu

    vc_list <- VarCorr(fit)$cond
    Gblocks <- lapply(vc_list, as.matrix)
    sigma2 <- sigma(fit)^2
    z <- numeric(length(e))

    for (g in seq_len(nlev)) {
        rows <- which(as.integer(group) == g)
        ng <- length(rows)
        V <- diag(sigma2, ng)
        for (i in seq_along(reStruc)) {
            bs <- reStruc[[i]]$blockSize
            cols <- Gp[i] + (g - 1) * bs + seq_len(bs)
            Zg <- as.matrix(Z[rows, cols, drop = FALSE])
            V <- V + Zg %*% Gblocks[[i]] %*% t(Zg)
        }
        R <- chol_with_jitter(V)
        z[rows] <- forwardsolve(t(R), e[rows])
    }

    list(
        term = vc_name,
        type = "marginal",
        group = names(flist)[1],
        fitted = mu,
        standardized = z
    )
}

##' Wrap Gaussian separable AR(1) models for DHARMa diagnostics
##'
##' \code{dharma_xar1} returns a lightweight wrapper around a fitted
##' \code{glmmTMB} model with a \code{homcsxar1} or \code{unxar1} random-effect
##' term. For Gaussian identity-link models, the wrapper lets
##' \code{DHARMa::simulateResiduals()} diagnose marginal normalized residuals,
##' i.e. fixed-effect residuals whitened by the fitted marginal covariance
##' implied by the model's random-effect terms.
##'
##' @param fit a fitted \code{glmmTMB} model containing a \code{homcsxar1} or
##' \code{unxar1} term
##' @param component unused; retained for compatibility
##' @param term optional \code{homcsxar1}/\code{unxar1} term name or numeric
##' index if the model contains more than one such term
##' @return an object suitable for \code{DHARMa::simulateResiduals()}
##' @details The Gaussian diagnostic is conditional on fitted parameter values
##' and uses the model-implied marginal covariance. Non-Gaussian models should
##' generally be checked with ordinary response-scale DHARMa simulations.
##' @examples
##' \dontrun{
##' res <- DHARMa::simulateResiduals(dharma_xar1(fit))
##' plot(res)
##' }
##' @export
dharma_xar1 <- function(fit, component = "cond", term = NULL) {
    if (component != "cond") {
        stop("dharma_xar1() currently supports only the conditional component.")
    }
    x <- get_xar1_marginal_residuals(fit, term = term)
    attr(fit, "dharma_xar1") <- x
    fit
}

xar1_dharma_data <- function(object) {
    attr(object, "dharma_xar1")
}

xar1_dharma_response <- function(object) {
    as.vector(xar1_dharma_data(object)$standardized)
}

xar1_dharma_fitted <- function(object) {
    as.vector(xar1_dharma_data(object)$fitted)
}

xar1_dharma_simulate <- function(object, nsim = 1, seed = NULL) {
    set_simcodes(object$obj, val = "random")
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
