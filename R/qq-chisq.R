 qq.chisq <-
  function(x, df=1, x.max,
    main="QQ plot", 
    sub=paste("Expected distribution: chi-squared (",df," df)", sep=""), 
    xlab="Expected", ylab="Observed",
    conc=c(0.025, 0.975), overdisp=FALSE, trim=0.5,  
    slope.one=FALSE, slope.lambda=FALSE, pvals=FALSE,
    thin=c(0.25,50), oor.pch=24, col.shade="gray", ...) {

    # Function to shade concentration band
    
    shade <- function(x1, y1, x2, y2, color=col.shade) {
      n <- length(x2)
      polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border=NA, col=color)
    }
    
    # Sort values and see how many out of range
    
    obsvd <- sort(x, na.last=NA)
    N <- length(obsvd)
    top <- obsvd[N]
    if (missing(x.max)) {
      Np <- N
    }
    else {
      Np <- sum(obsvd<=abs(x.max))
      if (Np<N || x.max<0) 
        top <- abs(x.max)
    }
    if(Np==0)
      stop("Nothing to plot")

    # Expected values
    
    if (df==2) {
      expctd <- 2*cumsum(1/(N:1))
    }
    else {
      expctd <- qchisq(p=(1:N)/(N+1), df=df)
    }

    # Concentration bands
    
    if (!is.null(conc)) {
      if(conc[1]>0) {
        e.low <- qchisq(p=qbeta(conc[1], 1:N, N:1), df=df)
      }
      else {
        e.low <- rep(0, N)
      }
      if (conc[2]<1) {
        e.high <- qchisq(p=qbeta(conc[2], 1:N, N:1), df=df)
      }
      else {
        e.high <- 1.1*rep(max(x),N)
      }
    }
    
    
    right <- expctd[N]
    par(mar=c(5,4,4,2)+0.1, las=1)
    if (pvals){
      mlp <- floor(-log10(pchisq(top, df=df, lower.tail=FALSE)))
      if (mlp>0) {
        gap <- ceiling(mlp/5)
        lp.vals <- seq(mlp, 1, -gap)
        chi2.vals <- qchisq(10^(-lp.vals), df=df, lower.tail=FALSE)
        par(mar=c(5,4,4,4)+0.1, las=1)
      }
      else 
        pvals <- FALSE
    }
    
    plot(c(0, right), c(0, top), type="n", xlab=xlab, ylab=ylab,
         main=main, sub=sub)

    # log p axis

    if (pvals) {
      nvals <- length(lp.vals)
      for (i in 1:nvals) 
        axis(side=4, at=chi2.vals[i],
             labels=substitute(10^{a}, list(a=-lp.vals[i])),  xaxt="n")
      mtext("P-value", side=4, line=3, las=0, padj=0)
    }
    
    
    # Thinning
    
    if (is.na(thin[1])) {
      show <- 1:Np
    }
    else if (length(thin)!=2 || thin[1]<0 || thin[1]>1 || thin[2]<1) {
      warning("invalid thin parameter; no thinning carried out")
      show <- 1:Np
    }
    else {
      space <- right*thin[1]/floor(thin[2])
      iat <- round((N+1)*pchisq(q=(1:floor(thin[2]))*space, df=df))
      if (max(iat)>thin[2]) 
        show <- unique(c(iat, (1+max(iat)):Np))
      else 
        show <- 1:Np
    }
    Nu <- floor(trim*N)
    if (Nu>0) 
      lambda <- mean(obsvd[1:Nu])/mean(expctd[1:Nu])
    if (!is.null(conc)) {
      if (Np<N)
        vert <- c(show, (Np+1):N)
      else
        vert <- show
      if (overdisp) 
        shade(expctd[vert], lambda*e.low[vert],
              expctd[vert], lambda*e.high[vert])
      else 
        shade(expctd[vert], e.low[vert], expctd[vert], e.high[vert])
    }
    points(expctd[show], obsvd[show], ...)
    # Overflow
    if (Np<N) {
      over <- (Np+1):N
      points(expctd[over], rep(top, N-Np), pch=oor.pch)
    }
    # Lines
    line.types <- c("solid", "dashed", "dotted")
    key <- NULL
    txt <- NULL
    if (slope.one) {
      key <- c(key, line.types[1])
      txt <- c(txt, "y = x")
      abline(a=0, b=1, lty=line.types[1])
    }
    if (slope.lambda && Nu>0) {
      key <- c(key, line.types[2])       
      txt <- c(txt, paste("y = ", format(lambda, digits=4), "x", sep=""))
      if (!is.null(conc)) {
        if (Np<N)
          vert <- c(show, (Np+1):N)
        else
          vert <- show
      }    
      abline(a=0, b=lambda, lty=line.types[2])
    }

    # Returned value
    
    if (!is.null(key)) 
       legend(0, top, legend=txt, lty=key)
    c(N=N, omitted=N-Np, lambda=lambda)

  }

    
