profilet.metaplus <- function(yi,sei,mods=NULL,justfit=FALSE,plotci=FALSE,slab=NULL) {
  
  isreg <- !is.null(mods)
  
  if (isreg) mods <- as.matrix(mods)
  
  hessll.profilet <- function(par,yi,sei,mods) {
    isreg <- !missing(mods)

     # browser()
    
    if (any(is.nan(par))) return(NA)
    
    muhat <- par[1]
    tau2 <- par[2]
    vinv <- par[3]
    if (isreg) xcoef <- matrix(par[4:length(par)],ncol=1)
    
    if (tau2 < 0.0) return(NA)
    
    f <- function(nu,oney,onesigma2,onemods,muhat,tau2,vinv,xcoef) {
      if ((tau2==0.0)|(vinv<0.0)) onel <- rep(0.0,length(nu))
      else {
        if (isreg) onel <- exp(-(oney-muhat-matrix(onemods,nrow=1) %*% xcoef-nu)^2/(2*onesigma2))*(1.0/sqrt(tau2))*dt(nu/sqrt(tau2),df=1.0/vinv)
        else onel <- exp(-(oney-muhat-nu)^2/(2*onesigma2))*(1.0/sqrt(tau2))*dt(nu/sqrt(tau2),df=1.0/vinv)
      }
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
    if (vinv==0.0){   
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
  
  ll.profilet <- function(par,yi,sei,mods) {
    negll <- hessll.profilet(par,yi,sei,mods)
    if (is.na(negll)) negll <- 1e100
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
  #browser()
  thenames <- c("muhat","tau2","vinv")
  if (isreg)  thenames <- c(thenames,dimnames(mods)[[2]])
  parnames(ll.profilet) <- thenames
  names(start.val) <- thenames
  names(lower.val) <- thenames
  
  # vinv=0.0 is a special case for some reason, maybe problem in optimisation
  normfit <- profilenorm.metaplus(yi=yi,sei=sei,mods=mods,justfit=TRUE,plotci=FALSE,slab=NULL)$fittedmodel
  if (isreg) start.null <- c(coef(normfit)[1:2],0.0,coef(normfit)[3:length(coef(normfit))])
  else start.null <- c(coef(normfit),0.0)
  names(start.null)[3] <- "vinv"
  #browser()
  if (isreg) maxfit <- mymle(ll.profilet,start=start.null,vecpar=TRUE,optimizer="user",
                                  data=list(yi=yi,sei=sei,mods=mods),
                                  skip.hessian=TRUE,
                                  control=list(eval.max=1000),
                                  lower=lower.val,optimfun=myoptim)
  else maxfit <- mymle(ll.profilet,start=start.null,vecpar=TRUE,optimizer="user",
                            data=list(yi=yi,sei=sei),
                            skip.hessian=TRUE,
                            control=list(eval.max=1000),
                            lower=lower.val,optimfun=myoptim)
  
  maxll <- logLik(maxfit)
  for (vinv in c(0.01,0.05,0.1,0.2,0.5,1)) { 
    start.val[3] <- vinv
    if (isreg) profilet.fit <- mymle(ll.profilet,start=start.val,vecpar=TRUE,optimizer="user",
                                                     data=list(yi=yi,sei=sei,mods=mods),
                                                     skip.hessian=TRUE,
                                                     control=list(eval.max=1000),
                                                     lower=lower.val,optimfun=myoptim)
    else profilet.fit <- mymle(ll.profilet,start=start.val,vecpar=TRUE,optimizer="user",
                                               data=list(yi=yi,sei=sei),
                                               skip.hessian=TRUE,
                                               control=list(eval.max=1000),
                                               lower=lower.val,optimfun=myoptim)
    if (logLik(profilet.fit) > maxll) {
      maxfit <- profilet.fit
      maxll <- logLik(profilet.fit)
    }
  }
  profilet.fit <- maxfit
# also tau2 and vinv are zero as not identifiable
  newstart.val <- start.val
  newstart.val[2:3] <- 0.0
  newlower.val <- lower.val
  newlower.val[2:3] <- 0.0
  if (isreg) profilet.fit2 <- mymle(ll.profilet,start=newstart.val,vecpar=TRUE,optimizer="user",
                                                     data=list(yi=yi,sei=sei,mods=mods),
                                                     skip.hessian=TRUE,
                                                     control=list(eval.max=1000),
                                                     lower=newlower.val,
                                                     fixed=list(tau2=0.0,vinv=0.0),optimfun=myoptim)
    else profilet.fit2 <- mymle(ll.profilet,start=newstart.val,vecpar=TRUE,optimizer="user",
                                               data=list(yi=yi,sei=sei),
                                               skip.hessian=TRUE,
                                               control=list(eval.max=1000),
                                                 lower=newlower.val,
                                               fixed=list(tau2=0.0,vinv=0.0),optimfun=myoptim)
  if(logLik(profilet.fit2) > logLik(profilet.fit)){
    profilet.fit@coef <- coef(profilet.fit2)
    profilet.fit@vcov <- vcov(profilet.fit2)
    profilet.fit@min <- profilet.fit2@min
    }
  results <- profilet.fit@coef
  if (!justfit) {
    if (isreg) thehessian <- hessian(hessll.profilet,results,yi=yi,sei=sei,mods=mods)
    else thehessian <- hessian(hessll.profilet,results,yi=yi,sei=sei)
    isproblem <- as.numeric(!(is.finite(diag(thehessian)) & (diag(thehessian)!=0.0)))
    isproblem2 <- isproblem*(1:length(results))
    noproblem2 <- (1-isproblem)*(1:length(results))
    if (all(isproblem2==0)) thehessian2 <- thehessian
    else thehessian2 <- thehessian[-isproblem2,-isproblem2]
    themyse <- suppressWarnings(sqrt(diag(ginv(thehessian2))))
    # expand back to original length
    myse <- rep(0,length(results))
    myse[noproblem2] <- themyse
    if (isreg) whichp <- c(1,4:(3+dim(mods)[2]))
    else whichp <- 1
    profilet.profile <- profilet.profile(profilet.fit,which=whichp,std.err=myse)
     profilet.ci <- confint(profilet.profile,method="uniroot")
    if (plotci) plot(profilet.profile)
    
    theci <- matrix(rep(NA,length(profilet.fit@coef)*2),ncol=2)
    theci[whichp,] <- profilet.ci
    results <- cbind(results,theci)
    pvalues <- rep(NA,length(start.val))
    for (iparm in whichp) {
      newstart.val <- coef(profilet.fit)
      newstart.val[iparm] <- 0.0
      newlower.val <- lower.val
      newlower.val[iparm] <- 0.0
      fixedparm <- names(start.val)[iparm]
      if (isreg) doprofile <- paste("profilet.fit0 <- mymle(ll.profilet,start=newstart.val,vecpar=TRUE,\n",
                                    "optimizer=\"user\",data=list(yi=yi,sei=sei,mods=mods),\n",
                                    "skip.hessian=TRUE,\n",
                                    "lower=newlower.val,fixed=list(",fixedparm,"=0.0),optimfun=myoptim)",sep="")
      else doprofile <- paste("profilet.fit0 <- mymle(ll.profilet,start=newstart.val,vecpar=TRUE,\n",
                              "optimizer=\"user\",data=list(yi=yi,sei=sei),\n",
                              "skip.hessian=TRUE,\n",
                              "lower=newlower.val,fixed=list(",fixedparm,"=0.0),optimfun=myoptim)",sep="")
       eval(parse(text=doprofile))
      pvalues[iparm] <- anova(profilet.fit,profilet.fit0)[2,5]
    }
    results <- cbind(results,pvalues)
    dimnames(results)[[2]] <- c("Est.","ci.lb","ci.ub","pvalue")
  }
  return(list(results=results,yi=yi,sei=sei,mods=mods,slab=slab,fittedmodel=profilet.fit,justfit=justfit,random="t-dist"))
}
