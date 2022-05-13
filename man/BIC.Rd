\name{BIC}
\alias{BIC.metaplus}
\title{BIC for metaplus object}
\description{
Returns BIC for a metaplus object.
}
\usage{
\method{BIC}{metaplus}(object, \ldots)
}
\arguments{
 \item{object}{metaplus object}
 \item{\dots}{ additional argument; currently none are used. }
}
\value{
BIC of fitted model
}
\author{Ken Beath}
\examples{
data(mag)
mag1 <- metaplus(yi, sei, plotci = TRUE, slab = study, cores = 1, data = mag)
print(BIC(mag1))
}
\keyword{methods}