\name{cdp}
\alias{cdp}
\docType{data}
\title{CDP meta-analysis data}
\description{
Data for the meta-analysis by Fioravanti and Yanagi (2005) of cytidinediphosphocholine (CDP-choline) for cognitive and behavioural disturbances associated with chronic cerebral disorders in the elderly.
}
\usage{data(cdp)}
\format{
  A data frame with 10 observations on the following 3 variables.
  \describe{
    \item{\code{study}}{study authors and date}
    \item{\code{yi}}{study effect estimate}
    \item{\code{sei}}{study standard error}
  }
}
\source{
  Fioravanti and Yanagi (2005)
  }
\references{
  Fioravanti, M., & Yanagi, M. (2005). Cytidinediphosphocholine (CDP choline) for cognitive and behavioural disturbances associated with chronic cerebral disorders in the elderly (Review). The Cochrane Database of Systematic Reviews. Retrieved from http://onlinelibrary.wiley.com/doi/10.1002/14651858.CD000269.pub2/pdf
}

\examples{
\dontrun{
data(cdp)
cdp1 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study)
}
}
\keyword{datasets}
