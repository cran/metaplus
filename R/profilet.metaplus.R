profilet.metaplus <- function(yi,sei,mods=NULL,justfit=FALSE,plotci=FALSE,slab=NULL) {
  
  isreg <- !is.null(mods)
  
  if (isreg) mods <- as.matrix(mods)
  
  ll.profilet <- function(par,yi,sei,mods) {
    isreg <- !missing(mods)
    
    if (any(is.nan(par))) return(NA)
    
    muhat <- par[1]
    tau2 <- par[2]
    vinv <- par[3]
    if (isreg) xcoef <- matrix(par[4:length(par)],ncol=1)
    
    if (tau2 < 0.0) return(NA)    
    if (vinv < 0.0) return(NA)
    
    
    f <- function(nu,oney,onesigma2,onemods,muhat,tau2,vinv,xcoef) {
         if (isreg) onel <- exp(-(oney-muhat-matrix(onemods,nrow=1) %*% xcoef-nu)^2/(2*onesigma2))*(1.0/sqrt(tau2))*dt(nu/sqrt(tau2),df=1.0/vinv)
        else onel <- exp(-(oney-muhat-nu)^2/(2*onesigma2))*(1.0/sqrt(tau2))*dt(nu/sqrt(tau2),df=1.0/vinv)
      return(onel)        
    }
    calconell <- function(x) {
      oney <- x[1]
      onesigma2 <- x[2]
      if (isreg) onemods <- x[3:length(x)]
      else onemods <- NULL
        theint <- tryCatch({
          if (isreg) integrate(f,-Inf,Inf,oney=oney,onesigma2=onesigma2,onemods=onemods,muhat=muhat,tau2=tau2,
                             vinv=vinv,xcoef=xcoef,subdivisions = 300L)$value
        else integrate(f,-Inf,Inf,oney=oney,onesigma2=onesigma2,onemods=NULL,muhat=muhat,tau2=tau2,vinv=vinv,
                        subdivisions = 300L)$value
      },
      error=function(e) {
        return(NA)
      })
      theint <- theint*(2*pi*onesigma2)^(-0.5)
      ll <- log(theint)
      return(ll)
    }
    if ((vinv==0.0) | (tau2==0.0)){   
      w <- 1.0/(tau2+sei^2)
      if (isreg) negll <- 0.5*sum(log(2*pi)+log(1/w)+w*(yi-muhat-as.vector(mods %*% xcoef))^2)
      else negll <- 0.5*sum(log(2*pi)+log(1/w)+w*(yi-muhat)^2)
    }
    else {
      if (isreg) negll <- -sum(apply(cbind(yi,sei^2,mods),1,calconell))
      else negll <- -sum(apply(cbind(yi,sei^2),1,calconell))
    }
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
     return(negll)
  }
  # obtain starting values
  if (isreg) {
    start.meta <- rma(yi=yi, sei=sei, mods=as.data.frame(mods), method="DL")
    start.val <- c(start.meta$b[1,1],start.meta$tau2,0.5,start.meta$b[2:dim(start.meta$b)[1],1])
    lower.val <- c(-Inf,0.0,0.0,rep(-Inf,dim(mods)[2]))
  } else {
    start.meta <- rma(yi=yi, sei=sei, method="DL")
    start.val <- c(start.meta$b[1,1],start.meta$tau2,0.5)
    lower.val <- c(-Inf,0.0,0.0)
  }
  thenames <- c("muhat","tau2","vinv")
  if (isreg)  thenames <- c(thenames,dimnames(mods)[[2]])
  parnames(ll.profilet) <- thenames
  names(start.val) <- thenames
  
#   # vinv=0.0 is a special case for some reason, maybe problem in nlminb?
#   maxfit <- profilenorm.metaplus(yi=yi,sei=sei,mods=mods,justfit=TRUE,plotci=FALSE,slab=NULL)$fittedmodel
#   if (isreg) maxfit@coef <- c(maxfit@coef[1:2],0.0,maxfit@coef[3:length(maxfit@coef)])
#   else maxfit@coef <- c(maxfit@coef[1:2],0.0)
#   names(maxfit@coef)[3] <- "vinv"
# #  print(logLik(maxfit))
#   maxll <- logLik(maxfit)
maxll <- -Inf
  for (vinv in c(0.0,0.01,0.05,0.1,0.2,0.5,1)) { 
    start.val[3] <- vinv
    if (isreg) profilet.fit <- suppressWarnings(mle2(ll.profilet,start=start.val,vecpar=TRUE,optimizer="nlminb",
                                                     data=list(yi=yi,sei=sei,mods=mods),
                                                     skip.hessian=TRUE,
                                                     control=list(eval.max=1000),
                                                     lower=lower.val))
    else profilet.fit <- suppressWarnings(mle2(ll.profilet,start=start.val,vecpar=TRUE,optimizer="nlminb",
                                               data=list(yi=yi,sei=sei),
                                               skip.hessian=TRUE,
                                               control=list(eval.max=1000),
                                               lower=lower.val))
    #print(logLik(profilet.fit))
    if (logLik(profilet.fit) > maxll) {
      maxfit <- profilet.fit
      maxll <- logLik(profilet.fit)
    }
  }
  profilet.fit <- maxfit
# also tau2 and vinv are zero as not identifiable
  newstart.val <- start.val[c(-2,-3)]
  newlower.val <- lower.val[c(-2,-3)]
  if (isreg) profilet.fit2 <- suppressWarnings(mle2(ll.profilet,start=newstart.val,vecpar=TRUE,optimizer="nlminb",
                                                     data=list(yi=yi,sei=sei,mods=mods),
                                                     skip.hessian=TRUE,
                                                     control=list(eval.max=1000),
                                                     lower=newlower.val,
                                                     fixed=list(tau2=0.0,vinv=0.0)))
    else profilet.fit2 <- suppressWarnings(mle2(ll.profilet,start=newstart.val,vecpar=TRUE,optimizer="nlminb",
                                               data=list(yi=yi,sei=sei),
                                               skip.hessian=TRUE,
                                               control=list(eval.max=1000),
                                                 lower=newlower.val,
                                               fixed=list(tau2=0.0,vinv=0.0)))
  if(logLik(profilet.fit2) >= logLik(profilet.fit)){
    profilet.fit@coef <- coef(profilet.fit2)
    profilet.fit@min <- profilet.fit2@min
    profilet.fit@details$convergence <- profilet.fit2@details$convergence    
    profilet.fit@details$message <- profilet.fit2@details$message
  } 
  if (profilet.fit@details$convergence!=0) {
# ignore singularity if due to small vinv
    if (!(((profilet.fit@details$message=="singular convergence (7)") & (profilet.fit@coef[3]<1e-6)) |
      (profilet.fit@details$message=="false convergence (8)") |
        (profilet.fit@details$message=="iteration limit reached without convergence (10)"))) {
      warning(paste("convergence failed: ",profilet.fit@details$message,sep=""))
    }
  } 
  results <- profilet.fit@coef
  if (!justfit) {
    if (isreg) thehessian <- hessian(ll.profilet,results,yi=yi,sei=sei,mods=mods)
    else thehessian <- hessian(ll.profilet,results,yi=yi,sei=sei)
    isproblem <- as.numeric(!(is.finite(diag(thehessian)) & (diag(thehessian)!=0.0)))
    isproblem2 <- isproblem*(1:length(results))
    noproblem2 <- (1-isproblem)*(1:length(results))
    if (all(isproblem2==0)) thehessian2 <- thehessian
    else thehessian2 <- thehessian[-isproblem2,-isproblem2]
    themyse <- sqrt(diag(solve(thehessian2)))
    # expand back to original length
    myse <- rep(0,length(results))
    myse[noproblem2] <- themyse
    if (isreg) whichp <- c(1,4:(3+dim(mods)[2]))
    else whichp <- 1
    profilet.profile <- suppressWarnings(profilet.profile(profilet.fit,which=whichp,std.err=myse))
     profilet.ci <- confint(profilet.profile,method="uniroot")
    if (plotci) plot(profilet.profile)
    
    theci <- matrix(rep(NA,length(profilet.fit@coef)*2),ncol=2)
    theci[whichp,] <- profilet.ci
    results <- cbind(results,theci)
    pvalues <- rep(NA,length(start.val))
    for (iparm in whichp) {
      newstart.val <- coef(profilet.fit)[-iparm]
      newlower.val <- lower.val[-iparm]
      fixedparm <- names(start.val)[iparm]
      if (isreg) doprofile <- paste("profilet.fit0 <- suppressWarnings(mle2(ll.profilet,start=newstart.val,vecpar=TRUE,\n",
                                    "optimizer=\"nlminb\",data=list(yi=yi,sei=sei,mods=mods),\n",
                                    "skip.hessian=TRUE,\n",
                                    "lower=newlower.val,fixed=list(",fixedparm,"=0.0)))",sep="")
      else doprofile <- paste("profilet.fit0 <- suppressWarnings(mle2(ll.profilet,start=newstart.val,vecpar=TRUE,\n",
                              "optimizer=\"nlminb\",data=list(yi=yi,sei=sei),\n",
                              "skip.hessian=TRUE,\n",
                              "lower=newlower.val,fixed=list(",fixedparm,"=0.0)))",sep="")
       eval(parse(text=doprofile))
      pvalues[iparm] <- anova(profilet.fit,profilet.fit0)[2,5]
    }
    results <- cbind(results,pvalues)
    dimnames(results)[[2]] <- c("Est.","ci.lb","ci.ub","pvalue")
  }
  return(list(results=results,yi=yi,sei=sei,mods=mods,slab=slab,fittedmodel=profilet.fit,justfit=justfit,random="t-dist"))
}
