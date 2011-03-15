convert.snpMatrix <- function(object) {
  cl <- class(object)
  pkg <- attr(cl, "package")
  newpkg<- "snpStats"
  old.classes <- c("snp.matrix", "X.snp.matrix", "single.snp.tests",
                   "single.snp.tests.score", "snp.tests.glm",
                   "snp.tests.glm.score", "snp.estimates.glm",
                   "imputation.rules")
  new.classes <- c("SnpMatrix", "XSnpMatrix", "SingleSnpTests",
                   "SingleSnpTestsScore", "GlmTests", "GlmTestsScore",
                   "GlmEstimates", "ImputationRules")
  icl <- match(cl[1], old.classes)
  if (is.na(icl)) {
    warning("unrecognized old snpMatrix class:", cl)
    return(object) ## return unchanged
  }
  else {
    if (icl==2) { ## X.snp.matrix
      object <- unclass(object)
      dip <- attr(object, "Female")
      object <- new("XSnpMatrix", diploid=dip)
      return(object)
    }
    else if (icl==5) { ## snp.tests.glm 
      object <-  new("GlmTests", snp.names=object@test.names,
                    var.names="Phenotype",
                    chisq=object@chisq, df=object@df, N=object@N)
      return(object)
    }
    else if (icl==6) { ## snp.tests.glm.score
      object <-  new("GlmTestsScore", snp.names=object@test.names,
                    var.names="Phenotype",                   
                    chisq=object@chisq, df=object@df, N=object@N,
                    score=object@score)
      return(object)
    }
    if (icl==8 && (is.null(pkg) || pkg!=newpkg)) 
      warning("Old snpMatrix imputation.rules will probably need to be regenerated")
    newcl <- new.classes[icl]
    attr(newcl, "package") <- newpkg
    class(object) <- newcl
    if (!isS4(object))
      object <- asS4(object)
    return(object)
  }
}

convert.snpMatrix.dir  <- function(dir=".", ext="RData") {
  files <- dir(dir)
  with.ext <- grep(paste(".", ext, "$", sep=""), files)
  if (length(with.ext)==0)
    cat("No files to be converted\n")
  else {
    toconvert <- files[with.ext]
    for (file in toconvert) {
      cat("Converting objects in file", file, "...")
      objects <- load(file)
      for (object in objects)
        assign(object, convert.snpMatrix(get(object)))
      cat("Saving")
      save(list=objects, file=file)
      cat("\n")
      remove(list=objects)
    }
  }
}




  
