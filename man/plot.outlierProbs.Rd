\name{plot.outlierProbs}
\alias{plot.outlierProbs}
\title{
Plot outlier probabilities.
}
\description{
Plots the outlier probability for each study, from an outlierProbs object.
}
\usage{
\method{plot}{outlierProbs}(x, ...)
}
\arguments{
 \item{x}{
outlierProbs object to be plotted
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
\donttest{
data(cdp)
cdp.mix <- metaplus(yi, sei, plotci = TRUE, slab = study, random = "mixture", cores = 1, data = cdp)
cdp.mix.outlierProbs <- outlierProbs(cdp.mix)
plot(cdp.mix.outlierProbs)
}
}
\keyword{method}
