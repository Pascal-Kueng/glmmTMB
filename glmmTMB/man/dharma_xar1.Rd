\name{dharma_xar1}
\alias{dharma_xar1}
\title{Wrap separable AR(1) latent innovations for DHARMa diagnostics}
\usage{
dharma_xar1(fit, component = "cond", term = NULL)
}
\arguments{
\item{fit}{a fitted \code{glmmTMB} model containing a \code{homcsxar1} or
\code{unxar1} term}

\item{component}{model component; currently defaults to the conditional
component}

\item{term}{optional \code{homcsxar1}/\code{unxar1} term name or numeric
index if the model contains more than one such term}
}
\value{
an object suitable for \code{DHARMa::simulateResiduals()}
}
\description{
\code{dharma_xar1} returns a lightweight wrapper around a fitted
\code{glmmTMB} model with a \code{homcsxar1} or \code{unxar1} random-effect
term. The wrapper lets \code{DHARMa::simulateResiduals()} diagnose the
standardized one-step latent AR(1) innovations rather than the observed
response. This is useful for Gaussian models fitted with
\code{dispformula = ~0}, where the lowest-level residual process is modeled
as a structured latent random effect.
}
\details{
The resulting diagnostics are based on fitted conditional modes of the latent
random effects. They are useful model checks, but should not be interpreted as
ordinary response-scale residual diagnostics.
}
\examples{
\dontrun{
res <- DHARMa::simulateResiduals(dharma_xar1(fit))
plot(res)
}
}
