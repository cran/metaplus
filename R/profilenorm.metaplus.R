profilenorm.metaplus <- function(yi,sei,mods=NULL,justfit=FALSE,plotci=FALSE,slab=NULL) {
  
  isreg <- !is.null(mods)
  
  ll.profilenorm <- function(par,yi,sei,mods) {
    isreg <- !missing(mods)
    muhat <- par[1]
    tau2 <- par[2]
    if (isreg) xcoef <- par[3:length(par)]
    w <- 1.0/(tau2+sei^2)
    if (isreg) negll <- 0.5*sum(log(2*pi)+log(1/w)+w*(yi-muhat-as.vector(as.matrix(mods) %*% xcoef))^2)
    else negll <- 0.5*sum(log(2*pi)+log(1/w)+w*(yi-muhat)^2)
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
    if (is.na(negll)) negll <- 1e100
    return(negll)
  }
    
  # obtain starting values
  
  if (isreg) {
    start.meta <- rma(yi=yi, sei=sei, mods=as.data.frame(mods), method="DL")
    start.vals <- c(start.meta$b[1,1],start.meta$tau2,start.meta$b[2:dim(start.meta$b)[1],1])
    lower.val <- c(-Inf,0.0,rep(-Inf,dim(mods)[2]))
  } else {
    start.meta <- rma(yi=yi, sei=sei, method="DL")
    start.vals <- c(start.meta$b[1,1],start.meta$tau2)
    lower.val <- c(-Inf,0.0)
  }
  if (isreg) names(start.vals) <- c("muhat","tau2",dimnames(mods)[[2]])
  else names(start.vals) <- c("muhat","tau2")
  
  parnames(ll.profilenorm) <- names(start.vals)
  
  names(lower.val) <- names(start.vals)
  
  start.vals <- start.vals+0.001
  
  if (isreg) profilenorm.fit <- mymle(ll.profilenorm,start=start.vals,vecpar=TRUE,
                                   optimizer="user",optimfun=myoptim,
                                   data=list(yi=yi,sei=sei,mods=as.matrix(mods)),
                                   skip.hessian=TRUE,
  #                                 control=list(eval.max=1000),
                                   lower=lower.val)
  else profilenorm.fit <- mymle(ll.profilenorm,start=start.vals,vecpar=TRUE,
                                optimizer="user",optimfun=myoptim,
                                data=list(yi=yi,sei=sei),
                             skip.hessian=TRUE,
 #                            control=list(eval.max=1000),
                             lower=lower.val)
  
  if (profilenorm.fit@details$convergence!=0) warning(paste("convergence failed: ",profilenorm.fit@details$message,sep="",))
  
  results <- profilenorm.fit@coef
 
  profilenorm.profiled <- NULL
  
  if (!justfit)  {
    notprofiled <- TRUE
    while (notprofiled) {
       if (isreg) thehessian <- hessian(ll.profilenorm,results,method.args=list(d=0.01),yi=yi,sei=sei,mods=as.matrix(mods))
      else thehessian <- hessian(ll.profilenorm,results,method.args=list(d=0.01),yi=yi,sei=sei)
      if (results[2] < 1.0e-6) {
        myse <- suppressWarnings(sqrt(diag(ginv(thehessian[-2,-2]))))
        if (length(myse)==1) myse <- c(myse,0.0)
        else myse <- c(myse[1],0.0,myse[2:length(myse)])
      } else myse <- suppressWarnings(sqrt(diag(ginv(thehessian))))
      
      if (isreg) whichp <- c(1,3:(2+dim(mods)[2]))
      else whichp <- 1
      profilenorm.profiled <- profile(profilenorm.fit,which=whichp,std.err=myse)
      if (inherits(profilenorm.profiled,"profile.mymle")) notprofiled <- FALSE
      else {
        thenames <- c("muhat","tau2")
        start.vals <- profilenorm.profiled@fullcoef
        if (isreg) {
          lower.val <- c(-Inf,0.0,rep(-Inf,dim(mods)[2]))
          thenames <- c(thenames,dimnames(mods)[[2]])
        } else {
          lower.val <- c(-Inf,0.0)
        }
        parnames(ll.profilenorm) <- thenames
        names(start.vals) <- thenames
        names(lower.val) <- thenames
        if (isreg) profilenorm.fit <- mymle(ll.profilenorm,start=start.vals,vecpar=TRUE,
                                           data=list(yi=yi,sei=sei,mods=mods),
                                           skip.hessian=TRUE,
                                           optimizer="user",optimfun=myoptim,
                                           control=list(eval.max=1000),
                                           lower=lower.val)
        else profilenorm.fit <- mymle(ll.profilenorm,start=start.vals,vecpar=TRUE,
                                     data=list(yi=yi,sei=sei),
                                      skip.hessian=TRUE,
                                     control=list(eval.max=1000),
                                     lower=lower.val,optimizer="user",optimfun=myoptim,)
        results <- profilenorm.fit@coef
      }
    }
    
    if (any(order(profilenorm.profiled@profile$muhat$z)!=(1:length(profilenorm.profiled@profile$muhat$z)))) 
      warning("Profile loglikelihood is not unimodal in region of estimate. Possibly incorrect confidence intervals.")
    profilenorm.ci <- confint(profilenorm.profiled,method="uniroot")
    if (plotci) {
      tryCatch(plot(profilenorm.profiled),
               error= function(e) {
                 print(paste("Error in CI plot: ",e))
               })
    }
    
    theci <- matrix(rep(NA,length(profilenorm.fit@coef)*2),ncol=2)
    theci[whichp,] <- profilenorm.ci
    results <- cbind(results,theci)
    
    # obtain p value for each coefficient
    pvalues <- rep(NA,length(start.vals))
    for (iparm in whichp) {
      fixedparm <- names(start.vals)[iparm]
      if (isreg) dostart <- paste("profilenorm.start <- makestart.profilenorm.metaplus(yi=yi,sei=sei,mods=as.matrix(mods),\n",
                                  "fixed=list(",fixedparm,"=0.0))",sep="")
      else dostart <- paste("profilenorm.start <- makestart.profilenorm.metaplus(yi=yi,sei=sei,\n",
                            "fixed=list(",fixedparm,"=0.0))",sep="")
      eval(parse(text=dostart))
      newstart.val <-  unlist(profilenorm.start)
      newstart.val <- newstart.val[names(start.vals)]
#      newstart.val <-  newstart.val[-iparm]
#      newlower.val <- lower.val[-iparm]
		newlower.val <- lower.val
      if (isreg) doprofile <- paste("profilenorm.fit0 <- mymle(ll.profilenorm,start=newstart.val,vecpar=TRUE,\n",
                                    "optimizer='user',optimfun=myoptim,data=list(yi=yi,sei=sei,mods=as.matrix(mods)),\n",
                                    "skip.hessian=TRUE,\n",
                                    "lower=newlower.val,fixed=list(",fixedparm,"=0.0))",sep="")
      else doprofile <- paste("profilenorm.fit0 <- mymle(ll.profilenorm,start=newstart.val,vecpar=TRUE,\n",
                              "optimizer='user',optimfun=myoptim,data=list(yi=yi,sei=sei),\n",
                              "skip.hessian=TRUE,\n",
                              "lower=newlower.val,fixed=list(",fixedparm,"=0.0))",sep="")
      eval(parse(text=doprofile))
      pvalues[iparm] <- anova(profilenorm.fit,profilenorm.fit0)[2,5]
    }
    results <- cbind(results,pvalues)
    dimnames(results)[[2]] <- c("Est.","95% ci.lb","95% ci.ub","pvalue")
  }
  return(list(results=results,yi=yi,sei=sei,mods=mods,slab=slab,justfit=justfit,fittedmodel=profilenorm.fit,profile=profilenorm.profiled,random="normal"))
}
