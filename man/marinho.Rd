\name{marinho}
\alias{marinho}
\docType{data}
\title{Marinho meta-analysis data}
\description{
Data for the meta-analysis by Marinho et al (2009) to determine the effectiveness of fluoride toothpastes on caries in children.
}
\usage{data(marinho)}
\format{
  A data frame with 70 observations on the following 11 variables.
  \describe{
    \item{\code{study}}{study authors and date}
    \item{\code{nfluor}}{number in fluoride group}
    \item{\code{meanfluor}}{mean effect in fluoride group}
    \item{\code{sdfluor}}{standard deviation of effect in fluoride group}
    \item{\code{nplacebo}}{number in placebo group}
    \item{\code{meanplacebo}}{mean effect in placebo group}
    \item{\code{sdplacebo}}{standard deviation of effect in placebo group}
    \item{\code{meaneffect}}{mean effect difference}
    \item{\code{seeffect}}{standard error of effect difference}
  }
}
\source{
  Marinho et al (2009)
  }
\references{
Marinho, V. C. C., Higgins, J. P. T., Logan, S., & Sheiham, A. (2009). Fluoride toothpastes for preventing dental caries in children and adolescents (Review). The Cochrane Database of Systematic Reviews. Retrieved from http://onlinelibrary.wiley.com/doi/10.1002/14651858.CD002278/pdf
}
\examples{
\dontrun{
data(marinho)
marinho1 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study)
}
}
\keyword{datasets}
