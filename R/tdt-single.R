#
# Much duplicated code between tdt.snp() and misinherits()
# Could/should be tidied up

tdt.snp <- function(ped, id, father, mother, affected,
                         data=sys.parent(), snp.data, rules=NULL,
                         snp.subset=NULL,
                         check.inheritance=TRUE, robust=FALSE, uncertain=FALSE,
                         score=FALSE) {
  if (!is.null(rules) || !check.inheritance) 
    robust <- TRUE
  mcall <- match.call()
  if (!is(snp.data, "SnpMatrix"))
    stop("snp.data argument must be of class SnpMatrix")
  if (missing(data)) { # ped data are in calling environment
    ped <- as.character(ped)
    nped <- length(ped)
    if (nped != nrow(snp.data))
      stop("length of `ped' argument incompatible with `snp.data'")
    id <- as.character(id)
    if (length(id)!=nped)
      stop("incompatible length for `id' and `ped' arguments")
    father <- as.character(father)
    if (length(father)!=nped)
      stop("incompatible length for `father' and `ped' arguments")
    mother <- as.character(mother)
    if (length(mother)!=nped)
      stop("incompatible length for `mother' and `ped' arguments")
    affected <- as.logical(affected)
    if (length(affected)!=nped)
      stop("incompatible length for `affected' and `ped' arguments")
               
    # Correspondence between ped data and snp data

    subject.names <- rownames(snp.data)
    in.snp <- 1:nped
    have.snps <- rep(TRUE, nped)
    
  } else { # ped data are in data dataframe
    data <- as.data.frame(data)
    nped <- nrow(data)
    subject.names <- rownames(data)

    if (missing(ped))
      ped <- as.character(data[,1])
    else
      ped <- as.character(eval(mcall$ped, envir=data))
    if (is.null(ped))
      stop("pedigree identifiers not found in data frame")

    if (missing(id))
      id <- as.character(data[,2])
    else
      id <- as.character(eval(mcall$id, envir=data))
    if(is.null(id))
      stop("subject identifiers not found in data frame")   

    if (missing(father))
      father <- as.character(data[,3])
    else
      father <- as.character(eval(mcall$father, envir=data))
    if(is.null(father))
      stop("father identifiers not found in data frame")   
    
    if (missing(mother))
      mother <- as.character(data[,4])
    else
      mother <- as.character(eval(mcall$mother, envir=data))
    if(is.null(mother))
      stop("mother identifiers not found in data frame")   
    
    if (missing(affected))
      affected <- (data[,6]==2)
    else
      affected <- as.logical(eval(mcall$affected, envir=data))
    if(is.null(affected))
      stop("disease status not found in data frame")   
  
    # Correspondence between ped data and snp data

    in.snp <- match(subject.names, rownames(snp.data))
    have.snps <- !is.na(in.snp)
  }

  # Treat subjects with "affected" missing as not affected
  
  affected[is.na(affected)] <- FALSE


  # Father and mother locations in ped file

  s.unique <- paste(ped, id, sep=":")
  if (any(duplicated(s.unique)))
    warning("Combination of pedigree ID and ID within pedigree does not generate unique IDs")
  f.unique <-  paste(ped, father, sep=":")
  fpos <- match(f.unique, s.unique)
  m.unique <-  paste(ped, mother, sep=":")
  mpos <- match(m.unique, s.unique)

  # Potentially complete trios
  
  trio <- have.snps & affected & (!is.na(fpos)) & (have.snps[fpos]) &
                                 (!is.na(mpos)) & (have.snps[mpos])
  ntrio <- sum(trio, na.rm=TRUE)
  if (ntrio==0) {
    cat("No potentially complete trios to analyse\n")
    return(NULL)
  }
  pd.snps <- in.snp[trio] # Proband rows in SnpMatrix
  fr.snps <- in.snp[fpos[trio]] # Fathers' rows in SnpMatrix
  mr.snps <- in.snp[mpos[trio]] # Mothers' rows in SnpMatrix
  
  clust <-   as.integer(factor(ped[trio]))
  cord <- order(clust)
  cat("Analysing", ntrio, "potentially complete trios in", max(clust),
      "different pedigrees\n")

  # Calculate scores and score variances

  scores <- .Call("score_tdt", pd.snps[cord], fr.snps[cord], mr.snps[cord],
                  clust[cord], snp.data, rules, snp.subset,
                  check.inheritance, robust, uncertain, PACKAGE="snpStats")
  chisq <- .Call("chisq_single", scores, PACKAGE="snpStats")
  if (is.null(rules)) {
    if (is.null(snp.subset))
      tested <- colnames(snp.data)
    else
      tested <- colnames(snp.data)[snp.subset]
  } else {
    if (is.null(snp.subset))
      tested <- names(rules)
    else
      tested <- names(rules)[snp.subset]
  }
  if (score)
    res <- new("SingleSnpTestsScore", snp.names=tested,
               chisq=chisq, N=scores$N, N.r2=scores$N.r2,
               U=scores$U, V=scores$V)
  else
    res <- new("SingleSnpTests", snp.names=tested, chisq=chisq,
               N=scores$N,  N.r2=scores$N.r2)
  res
}

misinherits <- function(ped, id, father, mother, data=sys.parent(), snp.data){

  non.mendel <- !as.logical(c(
     1,1,1,1, 1,1,1,0, 1,1,1,1, 1,0,1,1,
     1,1,1,0, 1,1,0,0, 1,1,1,0, 1,0,1,0,
     1,1,1,1, 1,1,1,0, 1,1,1,1, 1,0,1,1,
     1,0,1,1, 1,0,1,0, 1,0,1,1, 1,0,0,1,
     1,1,0,1, 1,1,0,0, 1,1,0,1, 1,0,0,1))

  mcall <- match.call()
  if (!is(snp.data, "SnpMatrix"))
    stop("snp.data argument must be of class SnpMatrix")
  X <- is(snp.data, "XSnpMatrix")
  if (missing(data)) { # ped data are in calling environment
    ped <- as.character(ped)
    nped <- length(ped)
    if (nped != nrow(snp.data))
      stop("length of `ped' argument incompatible with `snp.data'")
    id <- as.character(id)
    if (length(id)!=nped)
      stop("incompatible length for `id' and `ped' arguments")
    father <- as.character(father)
    if (length(father)!=nped)
      stop("incompatible length for `father' and `ped' arguments")
    mother <- as.character(mother)
    if (length(mother)!=nped)
      stop("incompatible length for `mother' and `ped' arguments")
    subject.names <- id
           
    # Correspondence between ped data and snp data

    subject.names <- rownames(snp.data)
    in.snp <- 1:nped
    have.snps <- rep(TRUE, nped)
           
  } else { # ped data are in data dataframe
    data <- as.data.frame(data)
    nped <- nrow(data)
    subject.names <- rownames(data)
    if (missing(ped))
      ped <- as.character(data[,1])
    else
      ped <- as.character(eval(mcall$ped, envir=data))
    
    if (missing(id))
      id <- as.character(data[,2])
    else
      id <- as.character(eval(mcall$id, envir=data))
    
    if (missing(father))
      father <- as.character(data[,3])
    else
      father <- as.character(eval(mcall$father, envir=data))
    
    if (missing(mother))
      mother <- as.character(data[,4])
    else
      mother <- as.character(eval(mcall$mother, envir=data))

    # Correspondence between ped data and snp data

    in.snp <- match(subject.names, rownames(snp.data))
    have.snps <- !is.na(in.snp)
  }

  # Father and mother locations in ped file

  s.unique <- paste(ped, id, sep=":")
  if (any(duplicated(s.unique)))
    warning("Combination of pedigree ID and ID within pedigree does not generate unique IDs")
  f.unique <-  paste(ped, father, sep=":")
  fpos <- match(f.unique, s.unique)
  m.unique <-  paste(ped, mother, sep=":")
  mpos <- match(m.unique, s.unique)

  # Potentially complete trios
  
  trio <- have.snps & (!is.na(fpos)) & (have.snps[fpos]) &
                                 (!is.na(mpos)) & (have.snps[mpos])
  ntrio <- sum(trio)
  if (ntrio==0) {
    cat("No potentially complete trios to analyse\n")
    return(NULL)
  }
  pd.snps <- in.snp[trio] # Proband rows in SnpMatrix
  fr.snps <- in.snp[fpos[trio]] # Fathers' rows in SnpMatrix
  mr.snps <- in.snp[mpos[trio]] # Mothers' rows in SnpMatrix

  fr.raw <- as.raw(snp.data[fr.snps,])
  if (X)
    fr.raw[!snp.data@diploid[pd.snps]] <- as.raw(4)
  mr.raw <- as.raw(snp.data[mr.snps,])
  pd.raw <- as.raw(snp.data[pd.snps,])
  ## Treat any uncertain genotypes as missing 
  fr.raw[fr.raw>3] <- as.raw(0)
  mr.raw[mr.raw>3] <- as.raw(0)
  pd.raw[pd.raw>3] <- as.raw(0)
  code <- 1 + as.numeric(
                rawShift(fr.raw, 4) |
                rawShift(mr.raw, 2) |
                pd.raw)
  error <- non.mendel[code]
  error[is.na(snp.data[pd.snps,])] <- NA
  error <- matrix(error, nrow=ntrio,
                  dimnames=list(subject.names[trio], colnames(snp.data)))
  erows <- apply(error, 1, any, na.rm=TRUE)
  ecols <- apply(error, 2, any, na.rm=TRUE)
  error[erows, ecols, drop=FALSE]
}

