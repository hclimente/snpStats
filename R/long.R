read.long <-
function(file, samples, snps,
         fields = c(snp = 1, sample = 2, genotype = 3, confidence = 4,
           allele.A=NA, allele.B=NA),
         split="\t| +", gcodes, no.call="", threshold=NULL,
         lex.order=FALSE, verbose=FALSE) {
  ## Number of ignored lines to list if verbose operation
  list.ignore <- 50
  ## Input file(s)
  if (length(file)>1) {
    if (missing(samples) || missing(snps))
      stop("When input is from multiple files, sample and snp arguments must be supplied")
    if (verbose)
      cat("Data to be read from multiple files", file[1], "...\n")
    con <- file(file[1], open="rt")
  } else {
    if (verbose)
      cat("Data to be read from the file", file, "\n")
    con <- file(file, open="rt")
  }
  ifile <- 1
  ## Sample and snp fields
  sampf <- fields["sample"]
  snpf <- fields["snp"]
  ## Size
  nsamp <- 0
  nsnp <- 0
  if (!missing(samples)) {
    if (is.character(samples)) {
      nsamp <- length(samples)
      if (verbose)
        cat("Names of", nsamp, "samples specified\n")
    }
    else if (is.numeric(samples) && length(samples)==1) {
      nsamp <- samples
      samples <- NULL
      if (verbose)
        cat("Number of samples to be read:", nsamp, "(names unspecified)\n")
    }
    else
      stop("illegal argument: samples")
  } else {
    samples <- NULL
  }
  if (!missing(snps)) {
    if (is.character(snps)) {
      nsnp <- length(snps)
      if (verbose)
        cat("Names of", nsnp, "snps specified\n")
    }    
    else if (is.numeric(snps) && length(snps)==1) {
      nsnp <- snps
      snps <- NULL
      if (verbose)
        cat("Number of snps to be read:", nsnp, "(names unspecified)\n")
    }
    else
      stop("illegal argument: snps")
  } else {
    snps <- NULL
  }

  ## Remaining argument checks
  
  cf <- fields["confidence"]
  if (is.null(threshold)) {
    if (!is.na(cf))
      stop("confidence field read but no thresholds set")
    if (verbose)
      cat("No confidence thresholds specified\n")
  } else {
    if (is.na(cf))
      stop("thresholds set but no confidence field read")
    low <- threshold[1]
    high <- threshold[2]
    if (verbose) {
      cat("Call rejected if confidence score ")
      if (!is.na(low))
        cat("is less than", low)
      if (length(threshold)>1)
        cat(" or ")
      if (!is.na(high))
        cat("is greater than", high)
      cat("\n")
    }       
  }
  af1 <- fields["allele.A"]
  af2 <- fields["allele.B"]
  gf <- fields["genotype"]
  if ((is.na(af1)+is.na(af2))==1)
    stop("only one allele field specified in fields argument")
  aread <- !is.na(af1) && !is.na(af2)
  gread <- !is.na(gf)
  if (aread) {
    if (gread)
      stop("both genotype and allele fields specified")
    if (verbose)
      cat("Alleles are read from two different fields\n")
    if (missing(gcodes)) {
      gcodes <- NA
    } else if (!is.na(gcodes)) {
      stop("gcodes argument set when data are read as alleles")
    }
  } else {
    if (!gread)
      stop("neither allele or genotype fields specified")
    if (verbose)
      cat("Genotype read as a single field ")
    if (missing(gcodes) || is.na(gcodes)) {
      ## 2-character allele coding
      if (verbose)
        cat("of two characters (which specify the alleles)\n")
      gcodes <- NULL
    } else if (length(gcodes)==3) {
      ## genotype coding
      gcodes <- as.character(gcodes)
      if (verbose)
        cat("matching one of", gcodes, "\n")
    } else if (!is.character(gcodes) || length(gcodes)!=1)
      stop("invalid value passed for gcodes argument")
    else if (verbose)
      cat("to be split into alleles by the regexp", gcodes, "\n")
  } 

  ## If missing size, read ahead

  if (nsamp==0 || nsnp==0) {
    if (verbose)
      cat("Initial scan of file\n")
    last.samp <- ""
    last.snp <- ""
    nline <- 0
    repeat {
      line <- readLines(con, n=1)
      if (length(line)==0 || nchar(line)==0) {
        if (nsamp==0) {
          nsamp <- length(samples)
          if (verbose)
            cat("Last sample:", last.samp, "\n")
        }
        if (nsnp==0) {
          nsnp <- length(snps)
          if (verbose)
            cat("Last snp:", last.snp, "\n")
        }
        break
      }
      nline <- nline+1
      fields <- strsplit(line, split)[[1]]
      this.samp <- fields[sampf]
      if (is.na(this.samp) || this.samp=="")
        stop("empty or NA sample field at line ", nline) 
      if (nsamp==0 && this.samp!=last.samp) {
        if (this.samp %in% samples) {
          if (verbose)
            cat("Last sample:", last.samp, "\n")
          nsamp = length(samples) ## Finished finding new sample names 
        } else {
          if (verbose && is.null(samples))
            cat("First sample:", this.samp, "\n")
          samples <- c(samples, this.samp)
          last.samp <- this.samp
        }
      }
      this.snp <- fields[snpf]
      if (is.na(this.snp) || this.snp=="")
        stop("empty or NA snp field at line ", nline)
      if (nsnp==0 && this.snp!=last.snp) {
        if (this.snp %in% snps) {
          if (verbose)
            cat("Last snp:", last.snp, "\n")
          nsnp = length(snps) ## Finished finding new snp names
        } else {
          if (verbose && is.null(snps))
            cat("First snp:", this.snp, "\n")
          snps <- c(snps, this.snp)
          last.snp <- this.snp
        }
      }
      if (nsamp>0 && nsnp>0)
        break;
    }
    if (nsamp==0 || nsnp==0)
      stop("Nothing read")
    seek(con, 0)  
  }
  
  ## Read data

  if (verbose)
    cat(nsamp, "x", nsnp, " matrix to be read\n", sep="")
  size <- nsamp*nsnp
  if (size>(2^31-1))
    stop("This is larger than the maximum size of a single object in R")
  genotypes <- matrix(as.raw(0), nrow=nsamp, ncol=nsnp)

  if (length(gcodes)==3)
    alleles <- NULL
  else
    alleles <- vector("list", nsnp)
  n.missing <- 0
  n.ignore <- 0
  n.reject <- 0
  n.found <- 0
  last.samp <- ""
  last.snp <- ""
  nline <- 0

  if (verbose) {
    cat("Reading genotypes from file\n ")
    for (i in 1:5) cat("       ", 2*i, "0%", sep="")
    cat("\n")
    for (i in 1:5) cat(".........|")
    cat("\n")
    bite <- size/50
    writ <- 0
  }
  
  while(n.found<size) {

    line <- readLines(con, n=1)
    if (length(line)==0 || nchar(line)==0) {
      if (ifile==length(file)) {
        break;
      } else {
        ifile <- ifile+1
        close(con)
        con <- file(file[ifile], open="rt")
        line <- readLines(con, n=1)
        if (length(line)==0 || nchar(line)==0)
          stop("Empty file: ", file[ifile])
      }
    }
    nline <- nline+1
    fields <- strsplit(line, split)[[1]]

    ## Sample
    
    this.samp <- fields[sampf]
    if (this.samp=="")
      stop("sample field empty at line ", nline)
    if (this.samp != last.samp) {
      i <- match(this.samp, samples)
      if (is.na(i)) {
        if (length(samples)<nsamp) {
          samples <- c(samples, this.samp)
          i <- length(samples)
        } else {
          n.ignore <- n.ignore+1
          next
        }
      }
      last.samp <- this.samp
    }

    ## Snp
    
    this.snp <- fields[snpf]
    if (this.snp=="")
      stop("snp field empty at line ", nline)
    if (this.snp != last.snp) {
      j <- match(this.snp, snps)
      if (is.na(j)) {
        if (length(snps)<nsnp) {
          snps <- c(snps, this.snp)
          j <- length(snps)
        } else {
          n.ignore <- n.ignore+1
          next
        }
      }
      last.snp <- this.snp
    }

    n.found <- n.found+1
    if (verbose && n.found>writ) {
      cat("-")
      writ <- writ+bite
    }
    ## Confidence filter

    if (!is.na(cf)) {
      fc <- fields[cf]
      if ((!is.na(low) && fc<low) || (!is.na(high) && fc>high)) {
        n.reject <- n.reject+1
        next
      }
    }

    ## Get genotype
    
    if (aread) {
      a12 <- c(fields[af1], fields[af2])
      if (a12[1]==no.call || a12[2]==no.call) {
        n.missing <- n.missing+1
        next
      }
    } else {
      fg <- fields[gf]
      if (fg==no.call) {
        n.missing <- n.missing+1
        next
      }
      if (is.null(gcodes)) {
        if (nchar(fg)!=2)
          stop("at line ", nline, ": ", fg,
               " (expecting a 2-character genotype field)")
        a12 <- c(substr(fg, 1, 1), substr(fg, 2, 2))
      } else if (length(gcodes)==3) {
        gtype <- match(fg, gcodes)
        genotypes[i,j] <- as.raw(gtype)
        next
      } else {
        a12 <- strsplit(fg, gcodes[1])
      }
    }
    aj <- alleles[[j]]
    if (length(aj)<2) {
      aj <- union(aj, a12)
      if (length(aj)>2)
        stop("on sample ", i, " snp ",  j, ":  not a SNP")
      alleles[[j]] <- aj
    }
    ia <- match(a12, aj)
    if (is.na(ia[1]) || is.na(ia[2]))
      stop("on sample ", i, " snp ",  j, ":  not a SNP")
    genotypes[i,j] <- as.raw(ia[1]+ia[2]-1)
  }
  if (verbose)
    cat("\n")
  close(con)
        
  if (n.ignore>0)
    warning(n.ignore, " lines were ignored")
  if (n.missing>0)
    warning(n.missing, " genotypes were not called")
  if (n.reject>0)
    warning(n.reject, " genotype calls did not pass confidence thresholds")
  
  ## Check that array is full
  
  msamp <- length(samples)
  msnp <- length(snps)

  if (msamp<nsamp || msnp<nsnp) {
    warning("Matrix was not filled; ", msamp, "x", nsnp, " matrix returned")
    genotypes <- genotypes[1:msamp, 1:msnp]
  }

  ## Make SnpMatrix object
  
  dimnames(genotypes) <- list(samples, snps)
  genotypes <- new("SnpMatrix", genotypes)

  ## Alleles frame, switching if necessary
  
  if (is.null(alleles)) {
    if (verbose)
      cat("Returning a SnpMatrix object\n")
    return(genotypes)
  } else {
    if (verbose)
      cat("Returning a list containing a SnpMatrix object and a dataframe listing the alleles\n")
    extract <- function(x, i) {
      if (is.null(x) || length(x)<i)
        return(as.character(NA))
      else
        return(x[i])
    }
    aA <- vapply(alleles, extract, character(1), 1)
    aB <- vapply(alleles, extract, character(1), 2)
    unseen <- is.na(aA) & is.na(aB)
    mono <- is.na(aB) & !unseen
    if (any(unseen))
      warning(sum(unseen), " snps had no acceptable calls")
    if (any(mono))
      warning(sum(mono), " snps were monomorphic")
    if (lex.order) {
      swa<- (!(is.na(aA)|is.na(aB)) & (aA>aB))
      if (any(swa)) {
        genotypes <- switch.alleles(genotypes, swa)
        aAsw <- aA[swa]
        aA[swa] <- aB[swa]
        aB[swa] <- aAsw
      }
    }
    alleles <- data.frame(allele.A=aA, allele.B=aB,
                          row.names=snps)
    return(list(genotypes=genotypes, alleles=alleles))
  }
}

