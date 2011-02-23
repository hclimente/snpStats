snp.imputation<- function(X, Y, pos.X, pos.Y, phase=FALSE, try=50,
                          stopping=c(0.95, 4, 0.05),
                          use.hap = c(1.0, 0.0),
                          em.cntrl=c(50, 0.01, 10, .01),
                          minA=5) {
  if (phase)
    stop("phase=TRUE option not yet implemented")
  if (!(is(X, "SnpMatrix") && is(Y, "SnpMatrix")))
    stop("X and Y arguments must be of class SnpMatrix or XSnpMatrix")
  if (class(X)[1]!=class(Y)[1])
    stop("X and Y arguments must be of the same class")
  if (nrow(X)!=nrow(Y) || length(pos.X)!=ncol(X) || length(pos.Y)!=ncol(Y))
    stop("Inconsistent argument dimensions")
  order.X <- order(pos.X)
  order.Y <- order(pos.Y)
  if (try>ncol(X))
    try <- ncol(X)
  .Call("snp_impute", X, Y, order.X, order.Y,
        as.double(pos.X[order.X]), as.double(pos.Y[order.Y]),
        as.logical(phase), as.integer(try),
        stopping, use.hap, em.cntrl, as.real(minA), PACKAGE="snpStats")
}

impute.snps <- function(rules, snps, subset=NULL, as.numeric=TRUE) {
  if (class(rules)!= "ImputationRules")
    stop("incorrect class for `rules' argument")
  if (class(snps)!= "SnpMatrix" && class(snps)!="XSnpMatrix")
    stop("incorrect class for `snps' argument")
  if (!is.null(subset)) {
    nr <- nrow(snps)
    if (is.logical(subset)) {
      if (length(subset)!=nr)
        stop("incompatible dimensions for snps and subset")
      subset <- (1:nr)[subset]
    }
    subset <- as.integer(subset);
    if (any((subset<1) | (subset>nr)))
      stop("row number out of range in subset argument")
  }
  .Call("impute_snps", rules, snps, subset, as.numeric, PACKAGE="snpStats")
}

# Returns names of rules (imputed snps) which feature excluded SNPs

filter.rules <- function(rules, snps.excluded, exclusions=TRUE) {
  if (class(rules)!= "ImputationRules")
    stop("incorrect class for `rules' argument")
  if (!is.character(snps.excluded))
    stop("excluded SNPs must be referenced by name")
  if.excl <- sapply(rules, FUN=function(rule, exc){
                             any(rule$snps %in% exc)
                           }, snps.excluded)
  if (exclusions)
    names(rules)[if.excl]
  else
    names(rules)[!if.excl]
}

# List extraction functions

imputation.maf <- function(rules) sapply(rules,
                    function(x) if (is.null(x$maf)) NA else x$maf)

imputation.r2 <- function(rules) sapply(rules,
                    function(x) if (is.null(x$r.squared)) NA else x$r.squared)

imputation.nsnp <- function(rules) sapply(rules,
                    function(x) if (is.null(x$snps)) NA else length(x$snps))
