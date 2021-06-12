.onAttach <-
  function (libname, pkgname) 
  {
    loadmsg <- "\nNote: The second parameter is the study standard error not the standard error squared as in the metafor package.\n"
    packageStartupMessage(loadmsg, domain = NULL, appendLF = TRUE)
  }

metaplus <- function(yi,sei,mods=NULL,random="normal",
      label=switch(random,"normal"="Random Normal","t-dist"="Random t-distribution", "mixture"="Random mixture"),
      plotci=FALSE,justfit=FALSE,slab=1:length(yi),
      useAGQ = FALSE,quadpoints=21,notrials=20, cores = max(detectCores() %/% 2, 1), 
      data) {
  if (!missing(useAGQ)) warning('useAGQ is deprecated as AGQ is no longer available')
  if (!missing(quadpoints)) warning('quadpoints is deprecated as AGQ is no longer available')
  if (!(random %in% c("normal","t-dist","mixture"))) stop("Unknown random effect type")
  if ((random=="mixture") & (notrials<10)) stop("Must be at least 10 sets of random starting values for mixture models.")
  if (cores<1) stop("Cores must be positive.")
  
  if (missing(data)) 
    data <- NULL
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  }
  else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  mf <- match.call()
  mf.yi <- mf[[match("yi", names(mf))]]
  mf.sei <- mf[[match("sei", names(mf))]]
  mf.slab <- mf[[match("slab", names(mf))]]
  mf.mods <- mf[[match("mods", names(mf))]]
  yi <- eval(mf.yi, data, enclos = sys.frame(sys.parent()))
  sei <- eval(mf.sei, data, enclos = sys.frame(sys.parent()))
  if (!is.null(mf.slab)) slab <- eval(mf.slab, data, enclos = sys.frame(sys.parent()))
  mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
  if (!is.null(mods)) {
    mods <- as.data.frame(mods)
    if (dim(mods)[2]==1) names(mods) <- deparse(mf.mods)
  }
  df <- switch(random,
                "normal"=length(yi)-1,
                "t-dist"=length(yi)-2,
                "mixture"=length(yi)-3)
  if (!is.null(mods)) df <- df-dim(mods)[2]
  if (df<=1) stop("Insufficient studies to fit model")
  if ((df<=3) & (!justfit)) warning("Very few studies. Solution may be unstable.")
  if (cores>1) loadNamespace("parallel")
  fit <- switch(random,
                "normal"=profilenorm.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab),
                "t-dist"=profilet.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab,quadpoints=quadpoints),
                "mixture"=profilemix.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab,notrials=notrials,cores=cores))
  fit$label <- label
  class(fit) <- "metaplus"
  return(fit)
}
  