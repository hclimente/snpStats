# Routines for structure analyses

xxt <- function(snps, strata=NULL, correct.for.missing=FALSE,
                lower.only=FALSE, uncertain = FALSE) {
  if (!is.null(strata) && !is.factor(strata))
    strata <- as.factor(strata)
  .Call("xxt", snps, strata, correct.for.missing, lower.only, uncertain,
        PACKAGE="snpStats")
}

ibsCount <- function(snps, uncertain=FALSE) {
  .Call("ibs_count", snps, uncertain, PACKAGE="snpStats")
}

ibsDist <- function(counts) {
  .Call("ibs_dist", counts, PACKAGE="snpStats")
}

snp.pre.multiply <- function(snps,  mat, frequency=NULL, uncertain=FALSE ) {
  .Call("snp_pre", snps, mat, frequency, uncertain, PACKAGE="snpStats")
}

snp.post.multiply <- function(snps,  mat, frequency=NULL, uncertain=FALSE) {
  .Call("snp_post", snps, mat, frequency, uncertain, PACKAGE="snpStats")
}

snp.cor <- function(x, y, uncertain=FALSE) {
  .Call("corsm", x, as.matrix(y), uncertain, PACKAGE="snpStats")
}
