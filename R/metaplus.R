metaplus <- function(yi,sei,mods=NULL,random="normal",
      label=switch(random,"normal"="Random Normal","t-dist"="Random t-distribution", "mixture"="Random mixture"),
      plotci=FALSE,justfit=FALSE,slab=1:length(yi)) {
  if (!(random %in% c("normal","t-dist","mixture"))) stop("Unknown random effect type")
  if (!is.null(mods)) mods <- as.data.frame(mods)
  df <- switch(random,
                "normal"=length(yi)-1,
                "t-dist"=length(yi)-2,
                "mixture"=length(yi)-3)
  if (!is.null(mods)) df <- df-dim(mods)[2]
  if (df<=1) stop("Insufficient studies to fit model")
  if ((df<=3) & (!justfit)) warning("Very few studies. Solution may be unstable.")
  fit <- switch(random,
                "normal"=profilenorm.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab),
                "t-dist"=profilet.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab),
                "mixture"=profilemix.metaplus(yi,sei,mods=mods,justfit=justfit,plotci=plotci,slab=slab))
  fit$label <- label
  class(fit) <- "metaplus"
  return(fit)
}
  