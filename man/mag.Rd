\name{mag}
\alias{mag}
\docType{data}
\title{Magnesium meta-analysis data}
\description{
Data for a meta-analysis of intravenous magnesium in acute myocardial infarction. An interesting question is whether the ISIS4 study is an outlier.
}
\usage{mag}
\format{
  A data frame with 16 observations on the following 3 variables.
  \describe{
    \item{\code{study}}{study author}
    \item{\code{yi}}{study effect estimate}
    \item{\code{sei}}{study standard error}
  }
}
\source{
  Sterne et al (2001)
  }
\references{
Sterne, J. A. C., Bradburn, M. J., & Egger, M. (2001). Meta-analysis in Stata. In M. Egger, G. D. Smith, & D. G. Altman (Eds.), Systematic Reviews in Health Care: Meta-Analysis in Context (pp. 347-369). BMJ Publishing Group.}

\examples{
\dontrun{
mag1 <- metaplus(yi,sei,plotci=TRUE,slab=study,data=mag)
}
}
\keyword{datasets}
