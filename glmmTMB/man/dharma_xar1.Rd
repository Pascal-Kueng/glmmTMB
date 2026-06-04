\name{dharma_xar1}
\alias{dharma_xar1}
\title{Wrap Gaussian separable AR(1) models for DHARMa diagnostics}
\usage{
dharma_xar1(fit, component = "cond", term = NULL)
}
\arguments{
\item{fit}{a fitted \code{glmmTMB} model containing a \code{homcsxar1} or
\code{unxar1} term}

\item{component}{unused; retained for compatibility}

\item{term}{optional \code{homcsxar1}/\code{unxar1} term name or numeric
index if the model contains more than one such term}
}
\value{
an object suitable for \code{DHARMa::simulateResiduals()}
}
\description{
\code{dharma_xar1} returns a lightweight wrapper around a fitted
\code{glmmTMB} model with a \code{homcsxar1} or \code{unxar1} random-effect
term. For Gaussian identity-link models, the wrapper lets
\code{DHARMa::simulateResiduals()} diagnose marginal normalized residuals,
i.e. fixed-effect residuals whitened by the fitted marginal covariance implied
by the model's random-effect terms.
}
\details{
The Gaussian diagnostic is conditional on fitted parameter values and uses the
model-implied marginal covariance. Non-Gaussian models should generally be
checked with ordinary response-scale DHARMa simulations.
}
\examples{
\dontrun{
res <- DHARMa::simulateResiduals(dharma_xar1(fit))
plot(res)
}
}
