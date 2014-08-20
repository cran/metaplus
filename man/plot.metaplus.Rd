\name{plot.metaplus}
\alias{plot.metaplus}
\title{
Produces forest plot for the studies together with the meta-analysis results.
}
\description{
Produces a forest plot for the studies in the meta-analysis and the result of the meta-analysis. Allows the inclusion of extra results of alternative meta-analyses, to allow, for example comparison between standard and robust methods of meta-analysis. Makes extensive use of the metafor package to produce the forest plot.}
\usage{
\method{plot}{metaplus}(x, ..., extrameta = NULL)
}
\arguments{
  \item{x}{
metaplus object to be plotted
}
 \item{\ldots}{
additional parameters to plot}
  \item{extrameta}{
Additional metaplus objects to be plotted. Only the summary statistic is plotted at the end of the plot. This is useful for comparing standard and robust methods, or meta-regressions for various values of the predictor.}
}
\value{
Plot
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\dontrun{
# produce all 3 models for the CDP data and plot them all
data(cdp)
cdp1 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study)
cdp2 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="t-dist")
cdp3 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="mixture")
plot(cdp1,extrameta=list(cdp2,cdp3))

# plot effect of exercise on depression at 4, 8 and 12 weeks
data(exercise)
exercise$duration4 <- exercise$duration-4
exercise$duration8 <- exercise$duration-8
exercise$duration12 <- exercise$duration-12

exercise.wk4 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration4"),
  drop=FALSE],label="Random Mixture (Week 4)",slab=exercise$study,random="mixture")
exercise.wk8 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration8"),
  drop=FALSE],label="Random Mixture (Week 8)",slab=exercise$study,random="mixture")
exercise.wk12 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration12"),
  drop=FALSE],label="Random Mixture (Week 12)",slab=exercise$study,random="mixture")

exercise.nodurn <- metaplus(exercise$smd,sqrt(exercise$varsmd),plotci=TRUE,
   label="Random Mixture (No Duration)",slab=exercise$study,random="mixture")

plot(exercise.nodurn,extrameta=list(exercise.wk4,exercise.wk8,exercise.wk12))
}
}
\keyword{method}
