row.summary <- function(object) {
   if (inherits(object, "SnpMatrix"))
     .Call("row_summary", object, PACKAGE="snpStats")
   else
     stop("not a SnpMatrix object")
}

col.summary <- function(object, rules=NULL, uncertain=TRUE) {
   if (inherits(object, "SnpMatrix")) {
     if (inherits(object, "XSnpMatrix")) 
       .Call("X_snp_summary", object, rules, uncertain, PACKAGE="snpStats")
     else  
       .Call("snp_summary", object, rules, uncertain, PACKAGE="snpStats")
   }
   else
     stop("not a SnpMatrix object")
 }

Fst <- function(snps, group) {
  .Call("Fst", snps, as.factor(group), PACKAGE="snpStats")
}

plotUncertainty <- function(snp, nlevels = 10,
                           color.palette = heat.colors(nlevels)) {
  codes <- factor(as.integer(snp), levels=1:253)
  freq <- table(codes)
  tot <- sum(freq, na.rm=TRUE)
  if (tot>0) {
    reord <- c(1, 253, 22, 2:21, 23:252)
    mass <- numeric(253)
    mass[reord] <- freq
    bym <- cut(mass, nlevels)
    h <- sqrt(3)/2
    d <- 1/(2*h)
    plot(c(0,21), c(0,23*h), type="n", xaxt="n", yaxt="n", xlab="", ylab="",
         bty="n", asp=1)
    lines(c(0, 21, 10.5, 0), c(0, 0, 21*h, 0))
    mtext("AA", side=1, line=0, at=0, adj=1, cex=1.5)
    mtext("BB", side=1, line=0, at=21, adj=0, cex=1.5)
    text(10.5, 22.3*h, "AB", cex=1.5, adj=0.5)
    ij <- 0
    for (i in 0:21) {
      yc <- i*h
      yll <- yc-d
      yl <- yc-0.5*d
      yh <- yc+0.5*d
      yhh <- yc+d
      xo <- i/2
      for (j in i:21) {
        ij <- ij+1
        xc <- j - xo
        xl <- xc-0.5   
        xr <- xc+0.5
        lev <- as.numeric(bym[ij])
        x <- c(xc, xl, xl, xc, xr, xr)
        y <- c(yll, yl, yh, yhh, yh, yl)
        mij <- mass[ij]
        if (mij>0) {
          polygon(x, y, col=color.palette[nlevels+1-lev], border=NA)
          text(xc, yc, as.character(mij), cex=0.6)
        }
      }
    }
  }
  else {
    "No data to plot"
  }
}
      
  
pp <- function(x, transpose=FALSE) {
  .Call("pp", as.raw(x), transpose, PACKAGE="snpStats")
}
