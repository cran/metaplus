testOutliers.profilemix.metaplus <- function(object,R=999,cores) {
  
  meta.fun <- function(data) {    
    
    isreg <- (dim(data)[2]>2)
     
    if (isreg) betadata.meta <-  metaplus(yi=data[,1], sei=data[,2], mods=data[,3:dim(data)[2],drop=FALSE],justfit=TRUE,plotci=FALSE,slab=NULL)
    else betadata.meta <-  metaplus(yi=data[,1], sei=data[,2],justfit=TRUE,plotci=FALSE,slab=NULL)
    if (isreg) betadata.meta2 <-  profilemix.metaplus(yi=data[,1], sei=data[,2], mods=data[,3:dim(data)[2],drop=FALSE],justfit=TRUE,plotci=FALSE,
                                                      cores=1,start.vals= c(betadata.meta$results[1],0.5*betadata.meta$results[2],max(0.1,1.5*betadata.meta$results[2]),log(0.05/(1-0.05)),betadata.meta$results[3:length(betadata.meta$results)]))
    else betadata.meta2 <-  profilemix.metaplus(yi=data[,1], sei=data[,2],justfit=TRUE,plotci=FALSE,
                                                cores=1,start.vals=c(betadata.meta$results[1],0.5*betadata.meta$results[2],max(0.1,1.5*betadata.meta$results[2]),log(0.05/(1-0.05))))
    if (abs((logLik(betadata.meta2$fittedmodel)-logLik(betadata.meta))/logLik(betadata.meta))<1.0e-4) sim.chisq <- 0
    else sim.chisq <- 2*(logLik(betadata.meta2$fittedmodel)-logLik(betadata.meta))
    if (sim.chisq < 0) {
      warning(paste("mixture chisq< 0","is",sim.chisq))
      sim.chisq <- 0
    }
    return(sim.chisq)
  }
  
  meta.rg <- function(data, mle) {
    isreg <- (dim(data)[2]>2)
    
    data[,2] <- sample(data[,2], length(data[,2]), replace = TRUE)
    if (isreg) mods <- data[,3:dim(data)[2],drop=FALSE]
    if (isreg) data[,1] <- rnorm(dim(data)[1],mean=mle$b[1,1],sd=sqrt(mle$tau2))+
      as.matrix(mods) %*%  mle$b[2:dim(mle$b)[1],1]+
      rnorm(dim(data)[1],mean=0,sd=data[,2])
    else data[,1] <- rnorm(dim(data)[1],mean=mle$b[1,1],sd=sqrt(mle$tau2))+
      rnorm(dim(data)[1],mean=0,sd=data[,2])
    return(data)
  }
  
  meta.ml <- rma(yi=object$yi, sei=object$sei, mods=object$mods, method="ML")

  if (cores>1) {
    if(.Platform$OS.type=="unix") parallel <- "multicore"
    else parallel <- "snow"
  } else parallel <- "no"
  
  meta.boot <- boot(data=cbind(object$yi,object$sei,object$mods), meta.fun, R = R, sim = "parametric",
                    ran.gen = meta.rg, mle=meta.ml, parallel=parallel,
                    ncpus = cores)
  pvalue <- (sum(meta.boot$t[,1] >= meta.boot$t0[1],na.rm=TRUE)+1)/(1+sum(!is.nan(meta.boot$t[,1])))
  return(list(pvalue=pvalue,observed=meta.boot$t0[1],sims=meta.boot$t[,1]))
}  