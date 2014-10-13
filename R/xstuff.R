.guessPloidy <- function(snps, known.diploid=NULL) {
  het <- row.summary(snps)$Heterozygosity
  N <- length(het)
  if (is.null(known.diploid)) {
    to.guess <- rep(TRUE, N)
    known.diploid <- logical(N)
  }
  else 
    to.guess <- is.na(known.diploid)
  dip <- ifelse(to.guess, het>0.1, known.diploid)
  repeat {
    diploid <- dip
    hp <- mean(het[diploid], na.rm=TRUE)
    dp <- mean(het[!diploid], na.rm=TRUE)
    if (is.na(dp) || is.na(hp))
      break
    hd <- (het-hp)^2/(hp*(1-hp))
    dd <- (het-dp)^2/((dp+0.001)*(1-dp))
    dip <- ifelse(to.guess, hd<dd, known.diploid)
    if (all(dip==diploid, na.rm=T))
      break
  }
  list(diploid=diploid, Heterozygosity=het)
}

.forceHom <- function(xsnps, diploid) {
  .Call("force_hom", xsnps, diploid, PACKAGE="snpStats")
}
