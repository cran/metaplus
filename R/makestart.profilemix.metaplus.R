makestart.profilemix.metaplus <- function(yi,sei,mods=NULL,fixed=NULL,notrials,cores) {
  
  fitonemcore <- function(yi,sei,outlierstarts,mods,fixed) {
 
    fitoneml2reg <- function(yi,sei,outliers,mods=NULL,fixed) {
      
      isreg <- !is.null(mods)
      
      rl2reg <- function(muhat,tau2,tau2out,lpoutlier,xcoef,yi,sei,mods,isreg,prop) {
        poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
        w <- 1.0/(tau2+sei^2)
        if (isreg) p1 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat-as.vector(as.matrix(mods) %*% xcoef))^2)
        else p1 <-  -0.5*(log(2*pi)-log(w)+w*(yi-muhat)^2)
        w <- 1.0/(tau2out+sei^2)
        if (isreg) p2 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat-as.vector(as.matrix(mods) %*% xcoef))^2)
        else p2 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat)^2)
        if (!missing(prop)) {
          ll <- prop*cbind(p1,p2)
          negll <- -sum(apply(l,1,sum))
        } else {
          l <- exp(cbind(log(1-poutlier)+p1,log(poutlier)+p2))
          negll <- -sum(log(apply(l,1,sum)))
        }
        if (is.nan(negll)) negll <- NA
        if (!is.finite(negll)) negll <- NA
        return(negll)
      }
      
      optimrl2reg <- function(p,lpoutlier,yi,sei,mods,isreg,prop,fixed) {
        p[names(fixed)] <- fixed
        return(rl2reg(p[1],p[2],p[3],lpoutlier,matrix(p[4:(3+ncoef)],ncol=1),yi,sei,mods,isreg,prop))
      }
      
      if (isreg) start.meta <- makestart.profilenorm.metaplus(yi=yi[outliers==0],sei=sei[outliers==0],mods=as.data.frame(mods[outliers==0,,drop=FALSE]),fixed=fixed)
      else start.meta <- makestart.profilenorm.metaplus(yi=yi[outliers==0],sei=sei[outliers==0],fixed=fixed)
      
      currtau2 <- start.meta$tau2
      poutlier <- sum(outliers)/length(outliers)
      if (poutlier<(0.5/length(outliers))) poutlier <- 0.5/length(outliers)
      currlpoutlier <- log(poutlier/(1-poutlier))
      currmuhat <- start.meta$muhat
      ncoef <- dim(mods)[2]
            
      if (isreg) {
        start.meta <- unlist(start.meta)
        currxcoef <- matrix(start.meta[3:(ncoef+2)],ncol=1)
      } else currxcoef <- NULL
 
       if (isreg) currtau2out <- sqrt(sum(outliers*(yi-currmuhat-as.vector(as.matrix(mods) %*% currxcoef))^2)/sum(outliers))
      else currtau2out <- sqrt(sum(outliers*(yi-currmuhat)^2)/sum(outliers))
      if (is.nan(currtau2out)) currtau2out <- 3*currtau2
      if (currtau2out==0.0) currtau2out <- 0.2
      
      # assemble into a vector to specify fixed
      if (length(fixed)>0) {
        current.vals <- c(currmuhat,currtau2,currtau2out,currlpoutlier,currxcoef)
        thenames <- c("muhat","tau2","tau2out","lpoutlier")
        if (isreg) thenames <- c(thenames,dimnames(mods)[[2]])
        names(current.vals) <- thenames
        current.vals[names(fixed)] <- fixed
        currmuhat <- current.vals[1]
        currtau2 <- current.vals[2]
        currtau2out <- current.vals[3]
        currlpoutlier <- current.vals[4]
        if (isreg) currxcoef <- matrix(current.vals[5:length(current.vals)],ncol=1)
      }
      
      currll <- -1.0e100
      nem <- 0
      repeat {
        nem <- nem+1
        # expectation step
        currpoutlier <- exp(currlpoutlier)/(1+exp(currlpoutlier))
        w <- 1.0/(currtau2+sei^2)
        if (isreg) ll1 <- -0.5*(log(2*pi)-log(w)+w*(yi-currmuhat-as.vector(as.matrix(mods) %*% currxcoef))^2)+log(1-currpoutlier)
        else ll1 <- -0.5*(log(2*pi)-log(w)+w*(yi-currmuhat)^2)+log(1-currpoutlier)
        w <- 1.0/(currtau2out+sei^2)
        if (isreg) ll2 <- -0.5*(log(2*pi)-log(w)+w*(yi-currmuhat-as.vector(as.matrix(mods) %*% currxcoef))^2)+log(currpoutlier)
        else ll2 <- -0.5*(log(2*pi)-log(w)+w*(yi-currmuhat)^2)+log(currpoutlier)
        
        l <- exp(cbind(ll1,ll2))
        prop <- l/apply(l,1,sum)
        # maximisation step
        # calculate outlier proportion
        poutlier <- sum(prop[,2])/dim(prop)[1]
        currlpoutlier <- log(poutlier/(1-poutlier))
        
        if (isreg) {
          startvals <- c(currmuhat,currtau2,currtau2out,currxcoef)
          names(startvals) <- c("muhat","tau2","tau2out",dimnames(mods)[[2]])
        } else {
          startvals <- c(currmuhat,currtau2,currtau2out)
          names(startvals) <- c("muhat","tau2","tau2out")
        }
        # ????? convert to Nelder_mead
        if (isreg) results.nlm <- nlminb(startvals,optimrl2reg,
                                         #               control=list(trace=6),
                                         lower = c(-Inf,0,0,rep(-Inf,ncoef)),
                                         lpoutlier=currlpoutlier,
                                         prop=prop,
                                         yi=yi,sei=sei,mods=mods,
                                         isreg=isreg,fixed=fixed)
        else results.nlm <- nlminb(startvals,optimrl2reg,
                                   # control=list(trace=6),
                                   lower = c(-Inf,0,0),
                                   lpoutlier=currlpoutlier,
                                   prop=prop,
                                   yi=yi,sei=sei,mods=NULL,
                                   isreg=isreg,fixed=fixed)
        currmuhat <- as.numeric(results.nlm$par)[1]
        currtau2 <- as.numeric(results.nlm$par)[2]
        currtau2out <- as.numeric(results.nlm$par)[3]
        if (isreg) currxcoef <- matrix(as.numeric(results.nlm$par)[4:(3+ncoef)],ncol=1)
        else currxcoef <- NULL
        
        lastll <- currll
        currll <- -rl2reg(currmuhat,currtau2,currtau2out,currlpoutlier,currxcoef,yi,sei,mods,isreg)
        if (abs((lastll-currll)/currll)<1.0e-6) break()
        if (nem >1000) break()
      }
      
      # make sure tau2out is greater than tau2
      if (currtau2 > currtau2out) {
        temp <- currtau2out
        currtau2out <- currtau2
        currtau2 <- temp
        currlpoutlier <- -currlpoutlier 
      }
      return(list(logLik=currll,params=c(currmuhat,currtau2,currtau2out,currlpoutlier,
                                            currxcoef)))
    }
    res <- fitoneml2reg(yi,sei,outlierstarts,mods,fixed)
    return(res)
  }
  
  thestatistic <- function(data, ...) {
        
    nmods <- dim(data)[2]-2

    if (nmods >0) {
      if (nmods==1) mods <- as.data.frame(cbind(data[, 3:dim(data)[2]]))
      else mods <- as.data.frame(data[, 3:dim(data)[2]])
      names(mods) <- names(data)[3:dim(data)[2]]
    }
   else mods <- NULL
    
    noutliers <- max(1,round(length(data$yi)*0.2))
    outliers <- matrix(sample(c(rep(1,noutliers),rep(0,length(data$yi)-noutliers)),length(yi)), nrow=1)
    
    onestat <- fitonemcore(data$yi,data$sei,outliers,mods,fixed)
 
    return(c(onestat$logLik, onestat$params))
  }
  isreg <- !is.null(mods)
  
  bootdata <- data.frame(yi=yi, sei=sei)
  
  if (isreg) bootdata <- cbind(bootdata, mods)
    
  infixed <- fixed
  fixed <- unlist(infixed)
  names(fixed) <- names(infixed)

  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  maxll <- -Inf 
 
    if (cores>1) {
     if(.Platform$OS.type=="unix")  parallel <- "multicore"
         else parallel <- "snow" }
        else parallel <- "no"
       theboot <- boot(bootdata, thestatistic, R=notrials, sim = "parametric",
         ran.gen = function(d, p) d,
         fixed=fixed,
         parallel = parallel,
         ncpus = cores)
    
    theboot <- theboot$t
    thefails <- is.na(theboot[,1])
    nfails <- sum(thefails)
    if (nfails>0) theboot <- theboot[!thefails,]
    if (nfails>0) warning(sprintf("Failed to obtain starting values for %i starting sets", nfails))
    
    theboot <- theboot[order(theboot[,1], decreasing=TRUE),]
    theboot <- theboot[1,]

    if (isreg) res <- list(muhat=theboot[2], tau2=theboot[3], tau2out=theboot[4],
                    lpoutlier=theboot[5], xcoef=matrix(theboot[6:length(theboot)], ncol=1))
    else res <- list(muhat=theboot[2], tau2=theboot[3], tau2out=theboot[4],
                lpoutlier=theboot[5])
    res <- list(logLik=theboot[1],params=res)
    res
}
