\name{summary}
\alias{summary}
\alias{summary.metaplus}
\title{
Summary of a metaplus object.}
\description{
Summarises the meta-analysis, giving the parameter estimates and goodness of fit statistics. For the robust methods this includes the values of the parameters relating to robustness.
}
\usage{
\method{summary}{metaplus}(object,...)
}
\arguments{
  \item{object}{
metaplus object to summarise.}
  \item{\dots}{
Additional parameters to summary.}
}
\value{
\item{results}{Matrix containing parameter estimates, confidence intervals and p values}
\item{fitstats}{List containing logLikelihood, AIC and BIC}
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
data(cdp)
cdp1 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study)
summary(cdp1)
}
\keyword{method}
