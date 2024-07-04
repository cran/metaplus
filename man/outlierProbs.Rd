\name{outlierProbs}
\alias{outlierProbs}
\alias{outlierProbs.metaplus}
\title{
Calculate outlier probabilities for each study.}
\description{
For the normal mixture random effect calculates the probability that each study is an outlier based on the posterior probability of it being an outlier.}
\usage{
\method{outlierProbs}{metaplus}(object)
}
\arguments{
 \item{object}{
A metaplus object with a mixture (robust) random effects distribution.
}
}
\value{
\item{outlier.prob}{Posterior probability that each study is an outlier}
\item{slab}{Labels corresponding to each study}
}
\details{
The outlier probabilities are obtained as the posterior probabilities of each study being an outlier based on the fitted mixture model.
}
\author{
Ken Beath <ken@kjbeath.id.au>
}
\examples{
\donttest{
data(mag)
mag3 <- metaplus(yi, sei, plotci = TRUE, slab = study, random = "mixture", cores = 1, data = mag)
mag3.outlierProbs <- outlierProbs(mag3)
plot(mag3.outlierProbs)
}
}
\keyword{method}
