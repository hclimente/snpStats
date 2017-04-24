test.allele.switch <- function(snps, snps2=NULL, split=NULL, prior.df=1) {
  if ((class(snps)!="SnpMatrix") && (class(snps)!="XSnpMatrix"))
    stop("argument `snps' is not a SnpMatrix or XSnpMatrix")
  if (is.null(snps2) && is.null(split))
    stop("either second SnpMatrix or split vector must be specified")
  if (is.null(snps2)) {
    if (length(split) != nrow(snps)) 
      stop("argument `split' not equal to number or rows in SnpMatrix")
    if (length(unique(split[!is.na(split)]))!=2)
      stop("argument `split' cannot be coerced into a 2-level factor")
    split <- as.factor(split)
  }
  else {
    if ((class(snps2)!="SnpMatrix") && (class(snps2)!="XSnpMatrix"))
      stop("argument `snps2' is not a SnpMatrix or XSnpMatrix")
    if (ncol(snps)!=ncol(snps2))
      stop("two SnpMatrix objects have different numbers of columns")
  }
  .Call("test_switch", snps, snps2, split, prior.df, PACKAGE="snpStats")
}

