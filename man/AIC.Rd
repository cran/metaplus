\name{AIC}
\alias{AIC.metaplus}
\title{AIC for metaplus object}
\description{
Returns AIC for a metaplus object.
}
\usage{
\method{AIC}{metaplus}(object, \ldots)
}
\arguments{
 \item{object}{metaplus object}
 \item{\dots}{ additional argument; currently none are used. }
}
\value{
AIC of fitted model
}
\author{Ken Beath}
\examples{
data(mag)
mag1 <- metaplus(yi, sei, plotci = TRUE, slab = study, cores = 1, data = mag)
print(AIC(mag1))
}
\keyword{methods}