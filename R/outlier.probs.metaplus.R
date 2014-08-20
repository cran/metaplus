outlier.probs <-
  ## Short form for generic function 
  function(object) UseMethod("outlier.probs")

outlier.probs.metaplus <- function(object) {
  if (!inherits(object, "metaplus"))
    stop("Use only with 'metaplus' objects.\n")
  outliers <- list(outlier.prob=object$outlier.prob,slab=object$slab)
  class(outliers) <- "outlier.probs"
  return(outliers)
}