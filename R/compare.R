sm.compare <- function(obj1, obj2, row.wise=TRUE, col.wise=TRUE) {
  if (!row.wise && !col.wise)
    stop("nothing to do")
  subjects <- intersect(rownames(obj1), rownames(obj2))
  if (length(subjects)==0)
    stop("objects have no rows in common")
  if (length(subjects)!=nrow(obj1) || length(subjects)!=nrow(obj2))
    warning("non-comformant rows; proceeding with ", length(subjects),
            " which are in common")
  snps <- intersect(colnames(obj1), colnames(obj2))
  if (length(snps)==0)
    stop("objects have no columns in common")
  if (length(snps)!=ncol(obj1) || length(snps)!=ncol(obj2))
    warning("non-comformant columns; proceeding with ", length(snps),
            " which are in common")
  ## Make into 3-D array, with object as last dimension
  x <- array(c(obj1[subjects,snps], obj2[subjects,snps]),
             dim=c(length(subjects), length(snps), 2))
  ## Compare function 
  comp <- function(mat) {
    c(sum(mat[,1]==mat[,2]),
      sum(mat[,1]!=mat[,2]),
      sum(mat[,1]==00 & mat[,2]==00),
      sum((mat[,1]==00 & mat[,2]!=00)|(mat[,1]!=00 & mat[,2]==00)),
      sum((mat[,1]==01 & mat[,2]==01)|(mat[,1]==03 & mat[,2]==03)),
      sum((mat[,1]==01 & mat[,2]==03)|(mat[,1]==03 & mat[,2]==01)),
      sum(mat[,1]==02 & mat[,2]==02),
      sum((mat[,1]==02 & mat[,2]!=02)|(mat[,1]!=02 & mat[,2]==02))
    )
  }
  snames <-  c("Agree", "Disagree", "NA.agree", "NA.disagree", "Hom.agree",
                    "Hom.switch", "Het.agree", "Het.Hom")
  ## Row compare
  if (row.wise) {
    rc <- t(apply(x, 1, comp))
    colnames(rc) <- snames
    rownames(rc) <- subjects
    if (!col.wise)
      return(rc)
  }
  ## Column compare
  if (col.wise) {
    cc <- t(apply(x, 2, comp))
    colnames(cc) <- snames
    rownames(cc) <- snps
    if (!row.wise)
      return(cc)
  }
  ## Both comparisons carried out
  list(row.wise=rc, col.wise=cc)
}
  

  
  
