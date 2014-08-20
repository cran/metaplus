\name{test.outliers}
\alias{test.outliers}
\alias{test.outliers.metaplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests for the presence of outliers.
}
\description{
For the $t$-distribution models this is a test that the degrees of freedom is infinitite, or equivalently that it's inverse is zero. For the mixture-normal model it is a test that the proportion of outliers is zero. As both tests involve a parameter on the boundary of the parameter space, asymptotic theory does not apply, so a parametric bootstrap is performed to determine the empirical distribution of the test statistic under the null hypothesis. The observed likelihood ratio statistic is then compared to this distribution to determine the p value.
}
\usage{
test.outliers(object, R = 999)
}
\arguments{
  \item{object}{
The meta-analysis for which the presence of outliers is to be tested.
}
  \item{R}{
Number of simulations (parametric bootstraps) used in testing the hypothesis. Initially, it may be useful to set this to a smaller value, to allow faster execution time.
}
}
\value{
pvalue=pvalue,observed=meta.boot$t0[1],sims
\item{pvalue}{p value obtained from parametric bootstrap}
\item{observed}{Observed value of the likelihood ratio test statistic}
\item{sims}{Simulated values of the test statistic under the null hypothesis}
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\dontrun{
data(cdp)
cdp3 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="mixture")
summary(test.outliers(cdp3))
}
}
\keyword{methods}
