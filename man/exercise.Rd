\name{exercise}
\alias{exercise}
\docType{data}
\title{Exercise meta-analysis data}
\description{
Lawlor & Hopker (2001) performed a meta-analysis of trials of exercise in
the management of depression, which was subsequently analysed using meta-
regression (Higgins & Thompson 2004) with duration of treatment as a covariate. There is a possible
outlier, the study by Reuter. While there are additional predictors, it seems excessive to use them given the small number of studies.
}
\usage{exercise}
\format{
 A data frame with 10 observations on the following 9 variables.
 \describe{
 \item{\code{study}}{study author}
 \item{\code{smd}}{study effect estimate}
 \item{\code{varsmd}}{study effect variance}
 \item{\code{sesmd}}{study effect standard error}
 \item{\code{abstract}}{study available as abstract only?}
 \item{\code{duration}}{length of study in weeks}
 \item{\code{itt}}{intention to treat analysis?}
 \item{\code{alloc}}{outcome assessor blinded}
 \item{\code{phd}}{phd thesis?}
 }
}
\source{
 Higgins and Thompson (2004)
 }
\references{
Higgins, J. P. T., & Thompson, S. G. (2004). Controlling the risk of spurious findings from meta-regression. Statistics in Medicine, 23(11), 166382. doi:10.1002/sim.1752

Lawlor, D. A., & Hopker, S. W. (2001). The effectiveness of exercise as an intervention in the management of depression: systematic review and meta-regression analysis of randomised controlled trials. BMJ, 322(31 March), 18.
}

\examples{
\donttest{
exercise1 <- metaplus(smd, sqrt(varsmd), mods = duration, slab = study,
    cores = 1, data = exercise)
exercise2 <- metaplus(smd, sqrt(varsmd), mods = cbind(duration, itt), slab = study,
    cores = 1, data = exercise)
}
}
\keyword{datasets}
