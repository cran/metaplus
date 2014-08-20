\name{metaplus-package}
\alias{metaplus-package}
\docType{package}
\title{
Fits random effects meta-analysis models including robust models
}
\description{
Allows fitting of random effects meta-analysis producing confidence intervals based on profile likelihood (Hardy and Thompson, 1996). Two methods of robust meta-analysis are included, based on either the $t$-distribution (Baker and Jackson (2008) and Lee and Thompson (2008)) or normal-mixture distribution (Beath, 2014). Tests can be performed for the need for a robust model, using a parametric bootstrap, and for the normal-mixture the identity of the outliers using the posterior probability. Plots are produced allowing a comparison between the results of each method. Where possible use has been made of the metafor package.

Future enhancements include the addition of binary meta-analysis using a likelihoiod approach with confidence intervals from profile likelihood.
}
\section{The metaplus function}{
This is the main function that allows fitting the models.metaplus objects may be plotted, using plot, tested for outliers using test.outliers. The results of tests.outliers may also be plotted.
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\dontrun{
data(marinho)
# fit standard normal random effects model
marinho1 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study)
summary(marinho1)
# fit t-distribution random effects model
marinho2 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study,
  random="t-dist")
summary(marinho2)
# test for presence of outliers
summary(test.outliers(marinho2))
# fit mixture normal random effects model
marinho3 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study,
  random="mixture")
summary(marinho3)
# test for presence of outliers
summary(test.outliers(marinho3))
# obtain posterior probability of each study being an outlier and plot
marinho3.outlier.probs <- outlier.probs(marinho3)
plot(marinho3.outlier.probs)
# produce forest plot including all 3 fitted models to compare
plot(marinho1,extrameta=list(marinho2,marinho3))
}
}
\references{
Baker, R., & Jackson, D. (2008). A new approach to outliers in meta-analysis. Health Care Management Science, 11(2), 121-131. doi:10.1007/s10729-007-9041-8

Beath, K. J. (2014). A finite mixture method for outlier detection and robustness in meta-analysis. Research Synthesis Methods, (in press). doi:10.1002/jrsm.1114

Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to meta-analysis with random effects. Statistics in Medicine, 15, 619-629.

Lee, K. J., & Thompson, S. G. (2008). Flexible parametric models for random effects distributions. Statistics in Medicine, 27, 418-434. doi:10.1002/sim
}
