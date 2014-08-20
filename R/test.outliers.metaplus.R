test.outliers <-
  ## Short form for generic function 
  function(object,R=999) UseMethod("test.outliers")

test.outliers.metaplus <- function(object,R=999) {
  if (!inherits(object, "metaplus"))
    stop("Use only with 'metaplus' objects.\n")

  if (object$justfit) stop("Cannot use with objects fitted with justfit=TRUE")
  
  if (object$random=="normal") stop("cant test for outliers with normal random effects model")
  test.outliers <- switch(object$random,
                "t-dist"=test.outliers.profilet.metaplus(object,R),
                "mixture"=test.outliers.profilemix.metaplus(object,R))
  class(test.outliers) <- "test.outliers"
  return(test.outliers)
}

summary.test.outliers <- function(object,...) {
  if (!inherits(object, "test.outliers"))
    stop("Use only with 'test.outliers' objects.\n")
  thesummary <- list(pvalue=object$pvalue,observed=object$observed)
  class(thesummary) <- "summary.test.outliers"
  thesummary
}

print.summary.test.outliers <- function(x,...) {
  if (!inherits(x, "summary.test.outliers"))
    stop("Use only with 'summary.test.outliers.metaplus' objects.\n")
  cat("Observed LRT statistic ")
  cat(sprintf("%1.1f",x$observed))
  cat(" p value ")
  thepvalue <- format.pval(x$pvalue, digits = 4, eps = 1.0e-4) 
  cat(sprintf("%s\n",thepvalue))
}
