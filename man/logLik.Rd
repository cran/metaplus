\name{logLik}
\alias{logLik.metaplus}
\title{log Likelikelihood for metaplus object}
\description{
Returns the log Likelihood for a metaplus object.
}
\usage{
\method{logLik}{metaplus}(object, \ldots)
}
\arguments{
 \item{object}{metaplus object}
 \item{\dots}{ additional argument; currently none are used. }
}
\value{
The loglikelihood of the fitted model.
}
\author{Ken Beath}
\examples{
data(mag)
mag1 <- metaplus(yi, sei, plotci = TRUE, slab = study, cores = 1, data = mag)
print(logLik(mag1))
}
\keyword{methods}