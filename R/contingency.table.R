# This routine takes a whole row of summary() into input
.contingency.table <- function(sum.case, sum.control) {
  OR <- sum.case$MAF / (1-sum.case$MAF) * (1-sum.control$MAF) /sum.control$MAF
  cat("OR: ", OR, "\n")
  
  case    <- c((1-sum.case$MAF), sum.case$MAF) * 2 * sum.case$Calls / sum.case$Call.rate
  control <- c((1-sum.control$MAF), sum.control$MAF) * 2 * sum.control$Calls / sum.control$Call.rate
  
  x <- cbind(control, case)
  # We'll just do fisher test, but it only does integers
  # print(fisher.test(x))
  x
}
