write.SnpMatrix <-

function(x, file, as.alleles=FALSE, append=FALSE, quote=TRUE, sep=" ", eol="\n",
         na="NA", row.names=TRUE, col.names=TRUE) {
  if (!is(x, "SnpMatrix"))
    stop("argument must be a SnpMatrix object")
  if (append && col.names)
    stop("col.names option is illegal in append mode")
  res <- .C("write_as_matrix", as.character(file),
     x@.Data, as.integer(nrow(x)), as.integer(ncol(x)),
     rownames(x), colnames(x), as.logical(as.alleles), as.logical(append),
     as.logical(quote), as.character(sep), as.character(eol),
     as.character(na), as.logical(row.names), as.logical(col.names),
     logical(1), PACKAGE="snpAssoc")
  error <- res[[15]]
  if (error==1)
    stop("Couldn't open output file")
  else
    c(nrow(x), ncol(x))
}
