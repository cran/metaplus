\name{plot.outlier.probs}
\alias{plot.outlier.probs}
\title{
Plot outlier probabilities.
}
\description{
Plots the outlier probability for each study, from an outlier.probs object.
}
\usage{
\method{plot}{outlier.probs}(x, ...)
}
\arguments{
  \item{x}{
outlier.probs object to be plotted
}
 \item{\ldots}{
additional parameters to plot}
}
\value{
Plot
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\dontrun{
data(cdp)
cdp.mix <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="mixture")
cdp.mix.outlier.probs <- outlier.probs(cdp.mix)
plot(cdp.mix.outlier.probs)
}
}
\keyword{method}
