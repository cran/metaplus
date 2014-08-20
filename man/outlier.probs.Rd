\name{outlier.probs}
\alias{outlier.probs}
\title{
Calculate outlier probabilities for each study.}
\description{
For the normal mixture random effect calculates the probability that each study is an outlier based on the posterior probability of it being an outlier.}
\usage{
outlier.probs(object)
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
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\dontrun{
data(cdp)
cdp3 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="mixture")
cdp3.outlier.probs <- outlier.probs(cdp3)
plot(cdp3.outlier.probs)
}
}
\keyword{method}
