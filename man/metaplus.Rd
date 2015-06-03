\name{metaplus}
\alias{metaplus}
\title{
Fits random effects meta-analysis models, using either a standard normal ditribution, a \eqn{t}-distribution or a mixture of normals for the random effect.
}
\description{
Allows fitting of random effects meta-analysis producing confidence intervals based on profile likelihood. Two methods of robust meta-analysis are included, based on either the \eqn{t}-distribution or normal-mixture distribution. 
}
\usage{
metaplus(yi,sei,mods=NULL,random="normal",
      label=switch(random,"normal"="Random Normal","t-dist"="Random t-distribution",
        "mixture"="Random mixture"),
      plotci=FALSE,justfit=FALSE,slab=1:length(yi))
}
\arguments{
  \item{yi}{
vector of observed effect size
}
  \item{sei}{
vector of observed standard errors
}
  \item{mods}{
data frame of covariates corresponding to each study
}
  \item{random}{
The type of random effects distribution. One of "normal","t-dist","mixture", for standard normal, \eqn{t}-distribution or mixture of normals respectively.
}
  \item{label}{
The label to be used for this model when plotting
}
  \item{plotci}{
Should profile be plotted for each confidence interval?
}
  \item{justfit}{
Should model only be fitted? If justfit=TRUE then profiling and likelihood ratio statistics are not calculated. Useful for when bootstrapping.
}
  \item{slab}{
Vector of character strings corresponding to each study.
}
}
\details{
As well as standard normal random effects meta-analysis, fitting of robust models incorporating either a \eqn{t}-distribution random effect (Lee and Thompson, 2008 or Baker and Jackson, 2008) or a mixture of normals (Beath, 2014) may be fitted. For all models the profile log-likelihood (Hardy and Thompson, 1996) is used to determine the confidence intervals for the parameters, with corresponding p values calculated using the likelihood ratio test, to give consistency between the confidence intervals and p values. The profile log-likelihood produces confidence intervals with better properties than REML or Der Simonian-Laird method (Brockwell and Gordon 2001, Jackson et al 2010)

For the mixture normal model, multiple starting values are used. This is achieved by starting with the no outlier model, then all possible single outlier models and choosing the model with the maximum likelihood. The process is repeated by adding all possible outliers to obtain the two outlier models and the maximum likelihood chosen. The process is repeated until the likelihood is not further maximised by adding outliers. This method gives a faster fitting time and is more reliable than using random assignments, as we have usually only a small number of outliers.
}
\value{
\item{results}{Matrix containing columns for estimate, lower 95\%, upper 95\% and p value, or if justfit=TRUE then only the parameter estimates}
\item{yi}{Effect sizes}
\item{sei}{Standard error of effect sizes}
\item{mods}{Modifiers for meta-regression}
\item{slab}{Study labels}
\item{justfit}{Value of justfit passed to the method.}
\item{fittedmodel}{final model returned by mle2}
\item{random}{Type of random effect}
\item{outlier.prob}{Outlier probabilities for robust mixture model only}
}
\references{
Baker, R., & Jackson, D. (2008). A new approach to outliers in meta-analysis. Health Care Management Science, 11(2), 121131.

Beath, K. J. (2014). A finite mixture method for outlier detection and robustness in meta-analysis. Research Synthesis Methods, 5, 285-293. doi:10.1002/jrsm.1114

Brockwell, S. E., & Gordon, I. R. (2001). A comparison of statistical methods for meta-analysis. Statistics in Medicine, 20, 825-840.

Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to meta-analysis with random effects. Statistics in Medicine, 15, 619-629.

Jackson, D., Bowden, J., & Baker, R. (2010). How does the DerSimonian and Laird procedure for random effects meta-analysis compare with its more efficient but harder to compute counterparts? Journal of Statistical Planning and Inference, 140, 961-970. doi:10.1016/j.jspi.2009.09.017

Lee, K. J., & Thompson, S. G. (2008). Flexible parametric models for random effects distributions. Statistics in Medicine, 27, 418-434. doi:10.1002/sim
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
data(cdp)
cdp.meta <- metaplus(cdp$yi,cdp$sei)
summary(cdp.meta)
}
\keyword{method}
