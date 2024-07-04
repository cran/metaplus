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
Ken Beath <ken@kjbeath.id.au>
}
\examples{
\donttest{
data(mag)
mag.mix <- metaplus(yi, sei, plotci = TRUE, slab = study, random = "mixture", cores = 1, data = mag)
mag.mix.outlierProbs <- outlierProbs(mag.mix)
plot(mag.mix.outlierProbs)
}
}
\keyword{method}
