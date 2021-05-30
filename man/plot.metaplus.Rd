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
\note{
On some systems the alignment of characters in the CI may be poor when saving the plots to pdf. A solution is to use the extrafont package and a fixed width font, for eample Courier New. Then if creating the plot using pdf() the parameter family = "Courier New" will use this font and fonts = "Courier New" will add it to the pdf.

An alternative is to use pdf.options() with the same parameters before the plot, and pdf.options(reset = TRUE) will produce the same effect.
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\donttest{
data(cdp)
# produce all 3 models for the CDP data and plot them all
cdp1 <- metaplus(yi, sei, plotci = TRUE, slab = study, cores = 1, data = cdp)
cdp2 <- metaplus(yi, sei, plotci = TRUE, slab = study, random = "t-dist", cores = 1, data = cdp)
cdp3 <- metaplus(yi, sei, plotci = TRUE, slab = study, random = "mixture", cores = 1, data = cdp)
plot(cdp1, extrameta = list(cdp2, cdp3))

# plot effect of exercise on depression at 4, 8 and 12 weeks
data(exercise)
exercise$duration4 <- exercise$duration-4
exercise$duration8 <- exercise$duration-8
exercise$duration12 <- exercise$duration-12

exercise.wk4 <- metaplus(smd, sqrt(varsmd), mods = duration4, 
 label = "Random Mixture (Week 4)", slab = study, random = "mixture", cores = 1, data = exercise)
exercise.wk8 <- metaplus(smd, sqrt(varsmd), mods = duration8, 
 label = "Random Mixture (Week 8)", slab = study, random = "mixture", cores = 1, data = exercise)
exercise.wk12 <- metaplus(smd, sqrt(varsmd), mods = duration12, 
 label = "Random Mixture (Week 12)", slab = study, random = "mixture", cores = 1, data = exercise)

exercise.nodurn <- metaplus(smd, sqrt(varsmd), plotci = TRUE, 
 label = "Random Mixture (No Duration)", slab = study, random = "mixture",
 cores = 1, data = exercise)

plot(exercise.nodurn, extrameta = list(exercise.wk4, exercise.wk8, exercise.wk12))
}
}
\keyword{method}
