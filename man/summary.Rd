\name{summary}
\alias{summary.metaplus}
\title{
Summary of a metaplus object.}
\description{
Summarises the meta-analysis, giving the parameter estimates and goodness of fit statistics. For the robust methods this includes the values of the parameters relating to robustness. In the summary the 95\% confidence intervals are calculated using the profile likelihood method. To produce corresponding p-values these are obtained from the likelihood ratio test statistic.
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
cdp1 <- metaplus(yi,sei,plotci=TRUE,slab=study,data=cdp)
summary(cdp1)
}
\keyword{method}
