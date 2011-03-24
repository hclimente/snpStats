ld <- function(x, y=NULL, depth=NULL, stats, symmetric=FALSE) {
  if (missing(stats))
    stop("LD stats must be specified")
  cstats <- c("LLR", "OR", "Q", "Covar", "D.prime", "R.squared", "R")%in%stats
  if (!any(cstats))
    stop("nothing to calculate")
  if (is.null(y)) {
    if (is.null(depth))
      stop("depth argument must be supplied")
    dmax <- ncol(x)-1
    if (depth>dmax) {
      depth <- as.integer(dmax)
      warning("depth too large; it has been reset to", dmax)
    }
  }
  .Call("ld", x, y, as.integer(depth), cstats, symmetric, PACKAGE="snpStats");
}
