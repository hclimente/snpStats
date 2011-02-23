.guessSex <- function(snps, known.female=NULL) {
  het <- row.summary(snps)$Heterozygosity
  N <- length(het)
  if (is.null(known.female)) {
    to.guess <- rep(TRUE, N)
    known.female <- logical(N)
  }
  else 
    to.guess <- is.na(known.female)
  fem <- ifelse(to.guess, het>0.1, known.female)
  repeat {
    female <- fem
    fm <- mean(het[female], na.rm=TRUE)
    mm <- mean(het[!female], na.rm=TRUE)
    if (is.na(fm) || is.na(mm))
      break
    fd <- (het-fm)^2/(fm*(1-fm))
    md <- (het-mm)^2/((mm+0.001)*(1-mm))
    fem <- ifelse(to.guess, fd<md, known.female)
    if (all(fem==female, na.rm=T))
      break
  }
  list(Female=female, Heterozygosity=het)
}

.forceHom <- function(xsnps, female) {
  .Call("force_hom", xsnps, female, PACKAGE="snpStats")
}
