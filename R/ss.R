#setOldClass(c("haplotype", "genotype"), test=TRUE)

#content should be raw - TODO: put restriction 
setClass("SnpMatrix", contains="matrix")

setClass("XSnpMatrix", representation("SnpMatrix", diploid="logical"),
         contains="SnpMatrix")
 
# constructor 
setMethod("initialize", 
          "SnpMatrix",
          function(.Object, ...) {
            if(is.null(.Object)) {
              stop("object is null to constructor");
            }
            .Object <- callNextMethod()
            # Please don't remove the tests - they are to do with promises and not redundant.
            # The accessor methods transparently passes to .Data, the assignment methods are different...
            if(!is.null(dim(.Object)) && (dim(.Object)[1] > 0) && (dim(.Object)[2] > 0)) {
              if (mode(.Object) != "raw") {
                cat("coercing object of mode ", mode(.Object), " to SnpMatrix\n")
                mode(.Object@.Data) <- "raw"
              }
              if (is.null(dimnames(.Object))) {
                cat("object has no names - using numeric order for row/column names\n")
                dimnames(.Object@.Data) <- c(list(as.character(1:(dim(.Object)[1]))), list(as.character(1:(dim(.Object)[2]))))
              }
            }
            .Object
          })

setMethod("initialize", 
          "XSnpMatrix", 
          function(.Object, ...) {
            .Object <- callNextMethod()
            if (length(.Object@diploid)==0){
              warning("diploid slot not supplied; it has been estimated from heterozygosity")
              guess <- .guessPloidy(.Object)
              .Object@diploid <- guess$diploid
              het <- guess$Heterozygosity
            }
            else {
              lend <- length(.Object@diploid)
              if (lend==1)
                .Object@diploid <- rep(.Object@diploid, nrow(.Object))
              else if (lend!=nrow(.Object))
                stop("diploid argument incorrect length")
              dpna <- is.na(.Object@diploid)
              if (any(dpna)) {
                warning("diploid argument contains NAs; these have been replaced by guesses based on heterozygosity")
                guess <- .guessPloidy(.Object, .Object@diploid)
                .Object@diploid <- guess$diploid
                het <- guess$Heterozygosity
              }
              else
                het <- row.summary(.Object)$Heterozygosity               
            }
            hetm <- !.Object@diploid & (!is.na(het)&(het>0))
            if (any(hetm)){
              warning(sum(hetm, na.rm=TRUE),
                      " heterozygous calls for haploid genotypes; these were set to NA")
              .Object@.Data <- .forceHom(.Object@.Data, .Object@diploid)
            }
            .Object
          })

setMethod("[", signature(x="SnpMatrix",i="ANY",j="ANY",drop="ANY"),
          function(x, i, j, drop=FALSE) {
            if (drop!=FALSE)
              stop("dimensions cannot be dropped from a SnpMatrix object")
            if (missing(i)) {
              if (missing(j))
                return(x)
              else
                x <-  x@.Data[,j,drop=FALSE]
            }
            else {
              if (missing(j))
                x <- x@.Data[i,,drop=FALSE]
              else 
                x <- x@.Data[i,j,drop=FALSE]
            }
            cl <- "SnpMatrix"
            attr(cl, "package") <- "snpStats"
            class(x) <- cl
            # setting an S4 class doesn't not automatically
            # set the object's internal tag; do it manually
            x <- asS4(x)
            x
          }
)

setMethod("[", signature(x="XSnpMatrix",i="ANY",j="ANY",drop="ANY"),
          function(x, i, j, drop=FALSE) {
            if (drop!=FALSE)
              stop("dimensions cannot be dropped from an XSnpMatrix object")
            diploid <- x@diploid
            if (missing(i)) {
              if (missing(j))
                return(x)
              else 
                x <-  x@.Data[,j, drop=FALSE]
            }
            else {
              if (is.character(i))
                i <- match(i, rownames(x))
              diploid <- diploid[i]
              if (missing(j))
                x <- x@.Data[i,,drop=FALSE]
              else 
                x <- x@.Data[i,j,drop=FALSE]
            }
            new("XSnpMatrix", x, diploid=diploid)
          })

# The sub-assignment methods

setMethod("[<-", signature(x="XSnpMatrix",i="ANY",j="ANY",value="XSnpMatrix"),
          function(x, i, j, value) {
            # The raw sub assignment should *just work*:
            if (!missing(i) & !missing(j)) {
              x@.Data[i,j] <- value
            } else if (missing(i) & missing(j)) {
              x@.Data[,] <- value
            } else if (missing(i)) {
              x@.Data[,j] <- value
            } else if (missing(j)) {
              x@.Data[i,] <- value
            }
            # All we want to do is to shuffle the rows
            # if a row index is passed
            if(!missing(i)) {
              slot(x, "diploid")[i] <- slot(value, "diploid")
            }
            x
          })

setAs("SnpMatrix", "numeric",
      function(from) {
        .Call("asnum", from, PACKAGE="snpStats")
      })

setAs("SnpMatrix", "character",
      function(from) {
        df <- dim(from)
        dnames <- dimnames(from)
        from <- 1+as.integer(from)
        to <- ifelse(from<5, c("NA", "A/A", "A/B", "B/B")[from], "Uncertain")
        dim(to) <- df
        dimnames(to) <- dnames
        to
      })

setAs("SnpMatrix", "XSnpMatrix",
      function(from) {
        new("XSnpMatrix", from)
      })

setAs("XSnpMatrix", "character",
      function(from) {
        df <- dim(from)
        ifr <- 1 + as.integer(from)
        offset <-  4*rep(from@diploid, ncol(from))
        to <- ifelse(ifr<5, c("NA", "A", "Error", "B",
                               "NA", "A/A", "A/B", "B/B")[ifr + offset],
                     "Uncertain")
        dim(to) <- df
        dimnames(to) <- dimnames(from)
        to
      })

setAs("matrix", "SnpMatrix",
      function(from) {
        if (is.numeric(from)) {
          valid  <- (from==0 | from==1 | from==2)
          if(!all(is.na(from) | valid)) 
            warning("values other than 0, 1 or 2 set to NA")
          to <- from+1
          to[is.na(from)|!valid] <- 0
          mode(to) <- "raw"
        }
        else if (is.raw(from)) {
          not.valid  <- (from>3)
          if (any(not.valid))
            warning("values other than 01, 02 or 03 set to NA (00)")
          to <- from
          to[not.valid] <- as.raw(0)
        }
        else 
          stop("Can only convert `numeric' or `raw' matrices") 
        new("SnpMatrix", to)
      })
           
# Summary of a SnpMatrix

setGeneric("summary")

setMethod("summary", "SnpMatrix",
   function(object) {
     list(rows = summary(row.summary(object)),
          cols = summary(col.summary(object, uncertain=TRUE)))
   })

setMethod("summary", "XSnpMatrix",
   function(object) {
     tab <- table(object@diploid)
     names(tab) <- ifelse(as.logical(names(tab)), "diploid", "haploid")
     list(sex = tab,
          rows = summary(row.summary(object)),
          cols = summary(col.summary(object, uncertain=TRUE)))
   })

# Imputation

setClass("ImputationRules", contains="list")

setMethod("summary", "ImputationRules",
   function(object) {
     info <- .Call("r2_impute", object, PACKAGE="snpStats")
     i2 <- info[,2]
     levs <- sort(unique(i2))
     labs <- paste(levs, "tags")
     size <- factor(i2, levels= levs, labels=labs)
     r2 <- info[,1]
     r2[r2>1] <- 1
     byr2 <- cut(r2, c((0:9)/10, 0.95, 0.99, 1.0), right=FALSE,
                 include.lowest=TRUE)
     table(byr2, size, dnn=c("R-squared", "SNPs used"), useNA="ifany")
   })

setMethod("[",
       signature(x="ImputationRules", i="ANY", j="missing", drop="missing"),
       function(x, i){
         if (is.character(i))
           i <- match(i, names(x), nomatch=0)
         ss <- x@.Data[i]
         names(ss) <- names(x)[i]
         res <- new("ImputationRules", ss)
         attr(res, "Max.predictors") <-
           max(sapply(res, function(rule) length(rule$snps)))
         res
       })


# Show and plot methods

setMethod("show", "ImputationRules",
   function(object) {
     to <- names(object)
     for (i in 1:length(object)) {
       if (is.null(object[[i]]$snps))
         cat(to[i], "~ No imputation available\n")
       else {
         if (is.null(object[[i]]$hap.probs)) 
           stop("Old imputation rule - not supported by this version")
         else 
           cat(to[i], " ~ ", paste(object[[i]]$snps, collapse="+"))
         cat(" (MAF = ", object[[i]]$maf, 
             ", R-squared = ", object[[i]]$r.squared,
             ")\n", sep="")
       }
     }
   })

setMethod("plot", signature(x="ImputationRules", y="missing"),
   function(x, y, ...) {
     mat <- summary(x)
     n <- nrow(mat)
     m <- ncol(mat)
     if (is.na(rownames(mat)[n])) {
       mat <- mat[-n,-m]
       n <- n-1
       m <- m-1
     }
     val <- barplot(t(mat[n:1,]), legend.text=TRUE, 
                    beside=F, col=heat.colors(m),
                    xlab=expression(r^2), ylab="Number of imputed SNPs",
                    names.arg=rep("", n), cex.lab=1.3, ...)
     mtext(rownames(mat)[n:1], at=val, side=1, las=2, line=0.5, cex=0.8)
    })

setMethod("show", "XSnpMatrix",
   function(object) {
     nr <- nrow(object)
     nc <- ncol(object) 
     dip <- object@diploid
     cat("An XSnpMatrix with ", nr, " rows (", sum(!dip), " haploid and ",
         sum(dip), " diploid), and ", nc, " columns\n", sep="")
     if (nr>1)
       cat("Row names: ", rownames(object)[1],"...", rownames(object)[nr],"\n")
     else
       cat("Row name: ", rownames(object)[1], "\n")
     if (nc>1)
       cat("Col names: ", colnames(object)[1],"...", colnames(object)[nc],"\n")
     else
       cat("Col name: ", colnames(object)[1],"\n")
       
   })

setMethod("show", "SnpMatrix",
   function(object) {
     nr <- nrow(object)
     nc <- ncol(object) 
     cat("A SnpMatrix with ", nr, "rows and ", nc,
         "columns\n")
     if (nr>1)
       cat("Row names: ", rownames(object)[1],"...", rownames(object)[nr],"\n")
     else
       cat("Row name: ", rownames(object)[1],"\n")
     if (nc>1)
       cat("Col names: ", colnames(object)[1],"...", colnames(object)[nc],"\n")
     else
       cat("Col name: ", colnames(object)[1],"\n")
   })


setMethod("is.na", "SnpMatrix", function(x){ x==0})

.rbind2 <- function(x,y){
  .External("snp_rbind",x, y, PACKAGE="snpStats")
}
snp.rbind <- function(...){
  .External("snp_rbind", ..., PACKAGE="snpStats")
}
.cbind2 <- function(x,y){
  .External("snp_cbind",x, y, PACKAGE="snpStats")
}
snp.cbind <- function(...){
  .External("snp_cbind", ..., PACKAGE="snpStats")
}

setMethod("rbind2", signature(x="SnpMatrix", y="SnpMatrix"), .rbind2)
setMethod("cbind2", signature(x="SnpMatrix", y="SnpMatrix"), .cbind2)



# Tests 

setClass("SingleSnpTests", 
         representation(snp.names="character", chisq="matrix",
                        N="integer", N.r2="numeric"))
setClass("SingleSnpTestsScore", 
         representation("SingleSnpTests", U="matrix", V="matrix"),
         contains="SingleSnpTests")
setClass("GlmTests",
         representation(snp.names="ANY", var.names="character",
                        chisq="numeric", df="integer",
                        N="integer"))
setClass("GlmTestsScore", representation("GlmTests", score="list"),
         contains="GlmTests")

              
setMethod("[",
          signature(x="SingleSnpTests", i="ANY",
                    j="missing", drop="missing"),
          function(x, i, j, drop) {
            if (is.character(i)) {
              i <- match(i, x@snp.names)
              if (any(is.na(i))) {
                warning(sum(is.na(i)),
                        " SNPs couldn't be found in SingleSnpTests object\n")
                i <- i[!is.na(i)]
              }
            } else if (is.logical(i)) {
              if (length(i)!=length(x@snp.names))
                stop("logical selection array has incorrect length")
              i <- (1:length(i))[i]
            } else if (is.numeric(i)) {
              if (min(i)<0 || max(i)>length(x@snp.names))
                stop("selection index out of range")
            } else {
              stop("illegal selection array")
            }
            if (length(x@N.r2)>0)
              N.r2 <- x@N.r2[i]
            else
              N.r2 <- numeric(0)
            new("SingleSnpTests",
              snp.names=x@snp.names[i], chisq=x@chisq[i,, drop=FALSE],
                N=x@N[i], N.r2=N.r2)
          })

setMethod("[",
          signature(x="SingleSnpTestsScore", i="ANY",
                    j="missing", drop="missing"),
          function(x, i, j, drop) {
            if (is.character(i)) {
              i <- match(i, x@snp.names)
              if (any(is.na(i))) {
                warning(sum(is.na(i)),
                        " SNPs couldn't be found in SingleSnpTests object\n")
                i <- i[!is.na(i)]
              }
            } else if (is.logical(i)) {
              if (length(i)!=length(x@snp.names))
                stop("logical selection array has incorrect length")
              i <- (1:length(i))[i]
            } else if (is.numeric(i)) {
              if (min(i)<0 || max(i)>length(x@snp.names))
                stop("selection index out of range")
            } else {
              stop("illegal selection array")
            }
            if (length(x@N.r2)>0)
              N.r2 <- x@N.r2[i]
            else
              N.r2 <- numeric(0)
            new("SingleSnpTestsScore",
              snp.names=x@snp.names[i], chisq=x@chisq[i,,drop=FALSE],
                N=x@N[i], N.r2=N.r2,
                U=x@U[i,,drop=FALSE], V=x@V[i,,drop=FALSE])
          })

setMethod("[",
          signature(x="GlmTests", i="ANY",
                    j="missing", drop="missing"),
           function(x, i, j, drop) {
             if (is.character(i))
               i <- match(i, names(x), nomatch=0)
             new("GlmTests",
                 snp.names = x@snp.names[i],
                 var.names = x@var.names,
                 chisq = x@chisq[i],
                 df = x@df[i],
                 N = x@N[i])
          })

setMethod("[",
          signature(x="GlmTestsScore", i="ANY",
                    j="missing", drop="missing"),
          function(x, i, j, drop) {
            if (is.character(i)) 
              i <- match(i, names(x), nomatch=0)
            new("GlmTestsScore",
                snp.names = x@snp.names[i],
                var.names = x@var.names,
                chisq = x@chisq[i],
                df = x@df[i],
                N = x@N[i],
                score = x@score[i])
          })
                   
setAs("SingleSnpTests", "data.frame",
      function(from) {
        if (length(from@N.r2)>0)
          data.frame(N=from@N, N.r2=from@N.r2,
                     Chi.squared=from@chisq,
                     P.1df=pchisq(from@chisq[,1], df=1, lower.tail=FALSE),
                     P.2df=pchisq(from@chisq[,2], df=2, lower.tail=FALSE),
                     row.names=from@snp.names)
        else
          data.frame(N=from@N,
                     Chi.squared=from@chisq,
                     P.1df=pchisq(from@chisq[,1], df=1, lower.tail=FALSE),
                     P.2df=pchisq(from@chisq[,2], df=2, lower.tail=FALSE),
                     row.names=from@snp.names)
      })

setMethod("summary", "SingleSnpTests",
          function(object) {
            summary(as(object, 'data.frame'))
          })

setMethod("show", "SingleSnpTests",
          function(object) {
            print(as(object, 'data.frame'))
          })

setAs("GlmTests", "data.frame",
      function(from) {
            chi2 <- from@chisq
            df <- from@df
            p <- pchisq(chi2, df=df, lower.tail=FALSE)
            N <- from@N
            data.frame(row.names=names(from),
                               Chi.squared=chi2, Df=df, p.value=p)
          })

setMethod("summary", "GlmTests",
          function(object) {
            summary(as(object, "data.frame"))
          })

setMethod("show", "GlmTests",
          function(object) {
            print(as(object, "data.frame"))
          })

# There are no standard generics for these, but the call to standardGeneric
# avoids a warning message 
  
setGeneric("p.value", function(x, df) standardGeneric("p.value"),
           useAsDefault=FALSE)

setGeneric("chi.squared", function(x, df) standardGeneric("chi.squared"),
           useAsDefault=FALSE)

setGeneric("deg.freedom", function(x) standardGeneric("deg.freedom"),
           useAsDefault=FALSE)

setGeneric("effect.sign", function(x, simplify) standardGeneric("effect.sign"),
           useAsDefault=FALSE)

setGeneric("sample.size", function(x) standardGeneric("sample.size"),
           useAsDefault=FALSE)

setGeneric("effective.sample.size",
           function(x) standardGeneric("effective.sample.size"),
           useAsDefault=FALSE)

setGeneric("pool2", function(x, y, score) standardGeneric("pool2"),
           useAsDefault=FALSE)

setGeneric("names")

setMethod("p.value", signature(x="SingleSnpTests", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            p <- pchisq(x@chisq[,df], df=df, lower.tail=FALSE)
            names(p) <- x@snp.names
            p
          })

setMethod("p.value", signature(x="GlmTests", df="missing"),
          function(x, df) {
            p <- pchisq(q=x@chisq, x@df, lower.tail=FALSE)
            p
          })


setMethod("chi.squared", signature(x="SingleSnpTests", df="numeric"),
          function(x, df) {
            if (df>2)
              stop("df must be 1 or 2")
            chi2 <- x@chisq[,df]
            names(chi2) <- x@snp.names
            chi2
          })

setMethod("chi.squared", signature(x="GlmTests", df="missing"),
          function(x, df) {
            chi2 <- x@chisq
            chi2
          })


setMethod("deg.freedom", signature(x="GlmTests"),
          function(x) {
            df <- x@df
            df
          })

setMethod("effect.sign", signature(x="GlmTests", simplify="logical"),
          function(x, simplify=TRUE) {
            res <- sapply(x@score, function(x) sign(x$U))
            res
         })

setMethod("effect.sign",
          signature(x="SingleSnpTestsScore", simplify="missing"),
          function(x, simplify) {
            res <- sign(x@U[,1])
            names(res) <- x@snp.names
            res
          })

setMethod("sample.size", signature(x="GlmTests"),
          function(x) {
            res <- x@N
            res
          })

setMethod("sample.size", signature(x="SingleSnpTests"),
          function(x) {
            res <- x@N
            names(res) <- x@snp.names
            res
          })

setMethod("effective.sample.size", signature(x="SingleSnpTests"),
          function(x) {
            if (length(x@N.r2)==0)
              res <- x@N
            else
              res <- x@N.r2
            names(res) <- x@snp.names
            res
          })

setMethod("names", signature(x="SingleSnpTests"), function(x) x@snp.names)

setMethod("names", "GlmTests",
          function(x) {
            objn <- x@snp.names
            if (is.list(objn)) {
              nobjn <- attr(objn, "names")
              if (!is.null(nobjn))
                nobjn
              else
                vapply(objn, function(x){
                  paste(x[1],x[length(x)], sep="...")}, "character")
            }
            else
              as.character(objn)
          })

setMethod("pool2",
          signature(x="SingleSnpTestsScore",y="SingleSnpTestsScore",
                    score="logical"),
          function(x, y, score) {
            all.snps <- union(x@snp.names, y@snp.names)
            can.pool <- intersect(x@snp.names, y@snp.names)
            x.only <- setdiff(x@snp.names, can.pool)
            y.only <- setdiff(y@snp.names, can.pool)
            if (length(can.pool)>0) {
              xsel <- match(can.pool, x@snp.names)
              ysel <- match(can.pool, y@snp.names)
              N <- x@N[xsel] + y@N[ysel]
              U <- x@U[xsel,,drop=FALSE] + y@U[ysel,,drop=FALSE]
              V <- x@V[xsel,,drop=FALSE] + y@V[ysel,,drop=FALSE]
              if (length(x@N.r2)>0) {
                if (length(y@N.r2)>0)
                  Nr2 <- x@N.r2[xsel] + y@N.r2[ysel]
                else
                  Nr2 <- x@N.r2[xsel] + y@N[ysel]
              } else {
                if (length(y@N.r2)>0)
                  Nr2 <- x@N[xsel]+  y@N.r2[ysel]
                else
                  Nr2 <- numeric(0)
              }
            } else {
              N <- NULL
              U <- NULL
              V <- NULL
              Nr2 <- numeric(0)
            }
            if (length(x.only)>0) {
              xsel <- match(x.only, x@snp.names)
              U <- rbind(U, x@U[xsel,,drop=FALSE])
              V <- rbind(V, x@V[xsel,,drop=FALSE])
              if (length(Nr2>0)) {
                if (length(x@N.r2)>0) 
                  Nr2 <- c(Nr2, x@N.r2[xsel])
                else
                  Nr2 <- c(Nr2, x@N[xsel])
              } else {
                if (length(x@N.r2)>0)
                  Nr2 <- c(N, x@N.r2[xsel])
              }
              N <- c(N, x@N[xsel])
            }
            if (length(y.only)>0) {
              ysel <- match(y.only, y@snp.names)
              U <- rbind(U, y@U[ysel,,drop=FALSE])
              V <- rbind(V, y@V[ysel,,drop=FALSE])
              if (length(Nr2>0)) {
                if (length(y@N.r2)>0) 
                  Nr2 <- c(Nr2, y@N.r2[ysel])
                else
                  Nr2 <- c(Nr2, y@N[ysel])
              } else {
                if (length(y@N.r2)>0)
                  Nr2 <- c(N, y@N.r2[ysel])
              }
              N <- c(N, y@N[ysel])
            }
            chisq <- .Call("chisq_single", list(U=U, V=V, N=N),
                           PACKAGE="snpStats")
            if (score)
              res <- new("SingleSnpTestsScore",
                         snp.names=c(can.pool, x.only, y.only),
                         chisq=chisq, N=N, U=U, V=V, N.r2=Nr2)
            else
              res <- new("SingleSnpTests",
                         snp.names=c(can.pool, x.only, y.only),
                         chisq=chisq, N=N, N.r2=Nr2)
            res
          })

# need to edit this

setMethod("pool2",
          signature(x="GlmTestsScore",y="GlmTestsScore",
                    score="logical"),
          function(x, y, score) {
            if (!is.null(x@var.names) && !is.null(y@var.names) &&
                (any(x@var.names!=y@var.names)))
              warning("tests to be pooled have non-matching var.names slots") 
            nm.x <- names(x)
            nm.y <- names(y)
            if (any(duplicated(nm.x))||any(duplicated(nm.y)))
              stop("At least one object has duplicated test names")
            to.pool <- intersect(nm.x, nm.y)
            if (length(to.pool)>0) {
              res <- .Call("pool2_glm",
                           x[to.pool],
                           y[to.pool], score, 
                           PACKAGE="snpStats")
 
            }
            else {
              res <- NULL
            }
            x.only <- setdiff(nm.x, to.pool)
            y.only <- setdiff(nm.y, to.pool)
            if (length(x.only)==0 && length(y.only)==0)
              return(res);
            ix <- match(x.only, nm.x, nomatch=0)
            iy <- match(y.only, nm.y, nomatch=0)
            if (is.null(res)) {
              res.chisq <- NULL
              res.df <- NULL
              res.N <- NULL
              res.score <- NULL
            }
            else {
              res.chisq <- res@chisq
              res.df <- res@df
              res.N <- res@N
              if (score)
                res.score <- res@score
            }
            if (score)
              res <- new("GlmTestsScore",
                         test.names=c(to.pool, x.only, y.only),
                         var.names=x@var.names,
                         chisq=c(res.chisq, x@chisq[ix], y@chisq[iy]),        
                         df=c(res.df, x@df[ix], y@df[iy]),        
                         N=c(res.N, x@N[ix], y@N[iy]),        
                         score=append(res.score,
                           append(x@score[ix], y@score[iy])))
            else
               res <- new("GlmTests",
                          test.names=c(to.pool, x.only, y.only),
                          var.names=x@var.names,
                          chisq=c(res.chisq, x@chisq[ix], y@chisq[iy]),        
                          df=c(res.df, x@df[ix], y@df[iy]),        
                          N=c(res.N, x@N[ix], y@N[iy]))       
            res
          })

# Switch allele methods 

setGeneric("switch.alleles",
           function(x, snps) standardGeneric("switch.alleles"),
           useAsDefault=FALSE)

setMethod("switch.alleles", signature(x="SnpMatrix", snps="ANY"),
          function(x, snps) {
            if (is.character(snps)) 
              snps <- match(snps, colnames(x))
            else if (is.logical(snps)) {
              if (length(snps)!=ncol(x))
                stop("logical snp selection vector is wrong length")
              snps <- (1:ncol(x))[snps]
            }
            else if (!is.integer(snps))
              stop("snp selection must be character, logical or integer")
            if (any(is.na(snps) | snps>ncol(x) | snps<1))
              stop("illegal snp selection")
            .Call("smat_switch", x, snps, PACKAGE="snpStats")
          })

setMethod("switch.alleles", signature(x="SingleSnpTestsScore", snps="ANY"),
          function(x, snps) {
            if (is.character(snps)) {
              snps <- match(snps, x@snp.names)
              if (any(is.na(snps)))
                warning(sum(is.na(snps)),
                        " SNP names were not found in tests object")
              snps <- snps[!is.na(snps)]
            } 
            ntest <- length(x@snp.names)
            if (is.logical(snps)) {
              if (length(snps)!=ntest)
                stop("incompatible arguments")
              if (sum(snps)==0)
                return(x)
            } else if (is.numeric(snps)) {
              if (length(snps)==0)
                return(x)
              if (max(snps)>ntest || min(snps)<1)
                stop("incompatible arguments")
            } else {
              stop("SNPs to be switched must be indicated by name, position, or by a logical vector")
            }
            res <- x
            res@U[snps,1] <- -x@U[snps,1]
            if (ncol(x@U)==3) {
              res@U[snps,2] <- -x@U[snps,2]
              res@U[snps,3] <- x@U[snps,3] - x@U[snps,2]
              res@V[snps,4] <- x@V[snps,4] - 2*x@V[snps,3] +
                x@V[snps,2]
              res@V[snps,3] <- x@V[snps,3] - x@V[snps,2]
            } else {
              res@U[snps,2] <- x@U[snps,2] - x@U[snps,1]
              res@V[snps,3] <- x@V[snps,3] - 2*x@V[snps,2] +
                x@V[snps,1]
              res@V[snps,2] <- x@V[snps,2] - x@V[snps,1]
            }
            res
          })

## This next will be slow but needed rarely

setMethod("switch.alleles", signature(x="GlmTestsScore",
                                      snps="character"),
          function(x, snps) {
            if (is.list(x@snp.names))
              switch <- lapply(x@snp.names, function(x, y) x%in%y, snps)
            else
              switch <- x@snp.names %in% snps
            N <- length(x@score)
            new.score <- vector("list", N)
            for (i in 1:N) {
              Ui <- x$U
              Vi <- x$V
              Si <- switch[i]
              lu <- length(Ui)
              ls <- length(Si)
              if (lu!=ls) {
                if (lu%%ls)
                  stop("error in GlmTest object")
                Si <- rep(Si, rep(lu/ls, ls))
              }
              new.score[[i]] <- list(U=ifelse(Si, -Ui, Ui), V=Vi)
            }
            new("GlmTestsScore", snp.names=x@snp.names, var.names=x@var.names,
                chisq=x@chisq, df=x@df, N=x@N,
                score=new.score)
          })


pool <- function(..., score=FALSE) {
  argl <- list(...)
  na <- length(argl)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na<2)
    stop("need at least two test sets to pool")
  if (na>2) {
    p2 <- pool2(..1, ..2, score=TRUE)
    r <- do.call(pool, c(p2, argl[3:na], score=score))
  }
  else
    r <- pool2(..1, ..2, score=score)
  r
}


## Parameter estimates

setClass("GlmEstimates", contains="list")

setMethod("[", signature("GlmEstimates", i="ANY", j="missing",
                         drop="missing"),
          function(x, i) {
            if (is.character(i))
              i <- match(i, names(x))
            res <- x@.Data[i, drop=FALSE]
            if (length(res)==0)
              return(NULL)
            names(res) <- names(x)[i]
            new("GlmEstimates", res)
           })

setMethod("show", "GlmEstimates",
          function(object) {
            len0 <- max(5, nchar(names(object)))
            len1 <- max(10,
                        sapply(object, function(x)
                               max(if(is.null(x)) 0 else nchar(x$Y.var))))
            len2 <- max(9,
                        sapply(object, function(x)
                               max(if(is.null(x)) 0 else nchar(names(x$beta)))))
            space <- 2
            nwidth <- 10
            di <- 5
            divide <- rep("-", len0+len1+len2+3*nwidth+5*space)
            cat("\n", format("Model", justify="right", width=len0),
                format("Y-variable", justify="right", width=len1+space),
                format("Parameter", justify="right", width=len2+space),
                format("Estimate", justify="right", width=nwidth+space),
                format("S.E.", justify="right", width=nwidth+space),
                format("z-value", justify="right", width=nwidth+space),"\n",
                divide, "\n", sep="")
            labs <- names(object)
            for (j in 1:length(object)) {
              labj <- labs[j]
              x <- object[[j]]
              if (is.null(x)) {
                cat(format(labj, justify="right", width=len1),
                    "- No estimates available\n") 
              }
              else {
                nyj <- x$Y.var
                nb <- length(x$beta)
                namep <- names(x$beta)
                ii <- 0
                for (i in 1:nb){
                  ii <- ii+i
                  bi <- x$beta[i]
                  si <- sqrt(x$Var.beta[ii])
                  zi <- bi/si
                  cat(format(labj, justify="right", width=len0), sep="")
                  cat(format(nyj, justify="right", width=len1+space), sep="")
                  labj <- ""
                  nyj <- ""
                  cat(format(namep[i], justify="right", width=len2+space),
                      format(bi, justify="right", width=nwidth+space,
                             digits=di),
                      format(si, justify="right", width=nwidth+space,
                             digits=di),
                      formatC(zi, format="f", width=nwidth+space, digits=3),
                      "\n", sep="")
                }
              }
              cat(divide, "\n", sep="")
            }
          }
          )

## Wald test

setAs("GlmEstimates", "GlmTests",
      function(from) {
        .Call("wald", from, package="snpStats")
      }
      )
        

## To do

## names method for snp and X.snp

