mvtests <- function(phenotype, sets, stratum, data=sys.parent(), snp.data,
                     rules=NULL, complete=TRUE, uncertain=FALSE, score=FALSE) {
  m <- match.call()
  tmiss <- missing(sets)
  if (tmiss)
    sets <- NULL
  smiss <- missing(stratum)
  if (smiss)
    stratum <- NULL
  if (!is(snp.data, "SnpMatrix"))
    stop("snp.data argument must be of class SnpMatrix")
  nr.snps = nrow(snp.data)
  if (missing(data)) { ## phenotype and stratum are in calling environment
    if (is.factor(phenotype)) 
      phenotype <- contrasts(phenotype)[as.numeric(phenotype),, drop=FALSE]
    if (!is.matrix(phenotype))
      stop("Phenotype must be a factor or a matrix")
    if (nrow(phenotype)!=nr.snps)
      stop("incompatible number of rows for phenotype vector")
    if (!smiss) {
      stratum <- as.integer(as.factor(stratum))
      if (length(stratum)!=nr.snps)
        stop("incompatible length for stratum vector")
      sorder <- order(stratum)
      sorder[is.na(stratum)] <- 0
    }
    else {
      sorder <- 1:nr.snps
    }
  }
  else { # phenotype, stratum, and subset  are in data dataframe
    nm.snps <- rownames(snp.data)
    nm.data <- rownames(data)
    nr.snps < nrow(snp.data)
    phenotype <- eval(m$phenotype, envir=data)
    if (is.factor(phenotype)) 
      phenotype <- contrasts(phenotype)[as.numeric(phenotype),]
    if (!is.matrix(phenotype))
      stop("Phenotype must be a factor or a matrix")
    # which points to rows in the data matrix
    which <- match(nm.snps, nm.data, nomatch=NA)
    phenotype <- phenotype[which,]
    if (nrow(phenotype)!=nr.snps)
      stop("incompatible number of rows for phenotype vector")
    if (!smiss) {
      stratum <- as.integer(as.factor(eval(m$stratum, envir=data)))[which]
      sorder <- order(stratum)
      sorder[is.na(stratum)] <- 0
    }
    else {
      sorder <- 1:nr.snps
    }
  }
  if (class(snp.data)=="XSnpMatrix") {
    stop("Not yet implemented for SNPs on X chromosome")
    if(any(is.na(snp.data@Female))) {
      warning("There are ", sum(is.na(snp.data@Female)),
              " subjects with NA for sex. These were ignored")
      sorder[is.na(snp.data@Female)] <- 0
    }
  }
    # Coerce sets argument to correct form #

  if(is.null(sets)) {
    sets <- as.integer(1:ncol(snp.data))
  } 
  else if (is.character(sets)) {
    sets <- .col.numbers(sets, colnames(snp.data))
  }
  else if (is.numeric(sets)) {
    sets <- as.integer(sets)
  }
  else if(is.list(sets)) {
    sets <- lapply(sets, .col.numbers, colnames(snp.data))
  } else {
    stop("illegal sets argument")
  }
  .Call("mvphen", phenotype, snp.data, rules, stratum, as.integer(sorder),
        sets, complete, uncertain, score, PACKAGE="snpAssoc")

}



                         
