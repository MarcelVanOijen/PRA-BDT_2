# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 2: Sampling-Based PRA

  library(mvtnorm)
  library(latex2exp)

  load( "data/l_xz.L.RData" )
  load( "data/l_xz.NL.RData" )
  
  PRA0 <- function( x, z, thr=0 ) {
    n    <- length(x)      ;  H       <- which(x < thr)  ;  n_H <-length(H)
    Ez_H <- mean( z[ H] )  ;  Ez_notH <- mean( z[-H] )
    pH   <- n_H / n        ;  V       <- Ez_notH - Ez_H  ;  R   <- pH * V
    return( c( pH=pH, V=V, R=R ) ) }

## Two Examples
### Example I: Linear Relationship

  set.seed(1)

  xz     <- l_xz.L[[1]] ; x <- xz[,1] ; z <- xz[,2]
  thr    <- -1
  PRA.L1 <- PRA0( x, z, thr )
  pH     <- PRA.L1["pH"] ; V <- PRA.L1["V"] ; R <- PRA.L1["R"]

  xz  <- l_xz.L[[1]] ; x <- xz[,1] ; z <- xz[,2]
  thr <- -1
  PRA0( x, z, thr )

  H  <- which(x < thr)
  Ez <- mean( z ) ; Ez_H <- mean( z[H] ) ; Ez_notH <- mean( z[-H] )

  par( mfrow=c(1,2) )
  plot  ( x, z, xlab="x", ylab="z", pch=".", xlim=c(-6,6), ylim=c(-4,4) )
  abline( v=thr, col="black", lty=2, lwd=2 )
  lines( x=c(thr,max(x)), y=rep(Ez_notH,2), col="red" , lty=1, lwd=1.5)
  lines( x=range(x)     , y=rep(Ez,2)     , col="blue", lty=1, lwd=1.5)
  lines( x=c(min(x),thr), y=rep(Ez_H,2)   , col="red" , lty=1, lwd=1.5)
  
  label.Ez_notH <- TeX(r"($\hat{E}[z|\neg H]$)")
  label.Ez      <- TeX(r"($\hat{E}[z]$)")
  label.Ez_H    <- TeX(r"($\hat{E}[z|H]$)")
  text( x=thr   , y=max(z) , adj=c(1.2, 0  ), "thr" )
  text( x=max(x), y=Ez_notH, adj=c(0.3,-0.1), label.Ez_notH, col="red"  )
  text( x=max(x), y=Ez     , adj=c(0.6, 1.3), label.Ez     , col="blue" )
  text( x=min(x), y=Ez_H   , adj=c(0.8, 1.2), label.Ez_H   , col="red"  )
  
  barplot( rbind( R, V-R ), beside=F, axisnames=F, asp=5 )
  mtext( c("V","R"), side=4, at=c(V,R), las=1 )

  thr.seq <- seq( -2, 2, length.out=101 )

  PRA.seq <- sapply( thr.seq, function(t){PRA0(x,z,t)} )
  pH.seq  <- PRA.seq["pH",] ; V.seq <- PRA.seq["V",] ; R.seq <- PRA.seq["R",]

  par( mfrow=c(1,3), mar=c(4,4,2,0) )

  plot  ( thr.seq, pH.seq, type="l",
          main="p[H]", xlab="thr", ylab="")
  plot  ( thr.seq, V.seq, col="black", type="l", ylim=c(0,max(V.seq,na.rm=T)),
          main="V, R", xlab="thr", ylab="" )
  points( thr.seq, R.seq, col="red"  , type="l" )
  legend( "bottomright", legend=c("V","R"), col=c("black","red"),
          lty=1, cex=0.75 )

### Example II: Nonlinear Relationship

  set.seed(1)

  xz      <- l_xz.NL[[1]] ; x <- xz[,1] ; z <- xz[,2]
  thr.seq <- seq( 0, 3, length.out=76 )
  PRA.seq <- sapply( thr.seq, function(t){PRA0(x,z,t)} )

  xz      <- l_xz.NL[[1]] ; x <- xz[,1] ; z <- xz[,2]
  thr.seq <- seq( 0, 3, length.out=76 )
  sapply( thr.seq, function(t){PRA0(x,z,t)} )

  par( mfrow=c(1,3), mar=c(4,4,2,0) )

  pH.seq <- PRA.seq["pH",] ; V.seq <- PRA.seq["V",] ; R.seq <- PRA.seq["R",]
  plot( x, z, xlab="x", ylab="z", asp=1, pch=".", main="Data" )
  plot  ( thr.seq, pH.seq, type="l",
          main="p[H]", xlab="thr", ylab="")
  plot  ( thr.seq, V.seq, col="black", type="l",
          ylim=c(0,max(V.seq,na.rm=T)),
          main="V, R", xlab="thr", ylab="" )
  points( thr.seq, R.seq, col="red"  , type="l" )
  legend( "bottomright", legend=c("V","R"), col=c("black","red"),
          lty=1, cex=1 )

## Uncertainty Quantification (UQ)
### UQ for $p[H]$

  pH_Be <- function( x, thr=0 ) {
    n  <- length(x) ; H    <- which(x < thr) ; n_H <- length(H)
    a  <- 1 + n_H   ; b    <- 1 + n - n_H
    pH <- a / (a+b) ; s_pH <- sqrt( a*b/(a+b+1) ) / (a+b)
    return( c( pH=pH, s_pH=s_pH ) )
  }

## PRA with UQ

  PRA <- function( x, z, thr=0 ) {
    n       <- length(x)     ; H         <- which(x < thr) ; n_H <- length(H)
    Ez_H    <- mean( z[ H] ) ; s_Ez_H    <- sqrt( var(z[ H]) /    n_H  )
    Ez_notH <- mean( z[-H] ) ; s_Ez_notH <- sqrt( var(z[-H]) / (n-n_H) )
    pH      <- n_H / n       ; V         <- Ez_notH - Ez_H ; R  <- pH * V
    s_pH    <- sqrt( pH*(1-pH) / n )
    s_V     <- sqrt( s_Ez_H^2 + s_Ez_notH^2 )
    s_R     <- sqrt( s_pH^2*s_V^2 + s_pH^2*V^2 + pH^2*s_V^2 )
    return( c(pH=pH,V=V,R=R,s_pH=s_pH,s_V=s_V,s_R=s_R) )
  }

  thr.seq <- seq(0,3,length.out=76)

  # All data of first nonlinear dataset
  xz      <- l_xz.NL[[1]]   ; x <- xz[,1] ; z <- xz[,2]
  PRA.seq <- sapply( thr.seq, function(t){PRA(x,z,t)} )
  pH.seq  <- PRA.seq["pH",] ; s_pH.seq <- PRA.seq["s_pH",]
  V.seq   <- PRA.seq["V" ,] ; s_V.seq  <- PRA.seq["s_V" ,]
  R.seq   <- PRA.seq["R" ,] ; s_R.seq  <- PRA.seq["s_R" ,]
  
  # Subsample of 10%
  n       <- length(x)      ; isub <- sample( 1:n, n/10 )
  xsub    <- x[isub]        ; zsub <- z[isub]
  PRA.sub <- sapply( thr.seq, function(t){PRA(xsub,zsub,t)} )
  pH.sub  <- PRA.sub["pH",] ; s_pH.sub <- PRA.sub["s_pH",]
  V.sub   <- PRA.sub[ "V",] ; s_V.sub  <- PRA.sub["s_V" ,]
  R.sub   <- PRA.sub[ "R",] ; s_R.sub  <- PRA.sub["s_R" ,]

  par( mfrow=c(2,3), mar=c(4,4,2,0) )

  plot  ( x, z, xlab="", ylab="z", asp=1, pch=".", main="Data" )
  plot  ( thr.seq, pH.seq, type="l", main="p[H]", xlab="", ylab="")
  points( thr.seq, pH.seq+2*s_pH.seq, type="l", lty=2 )
  points( thr.seq, pH.seq-2*s_pH.seq, type="l", lty=2 )
  plot  ( thr.seq, V.seq, col="black", type="l",
          ylim=c(0,max(V.seq,na.rm=T)), main="V, R", xlab="", ylab="" )
  points( thr.seq, V.seq+2*s_V.seq, type="l", lty=2 )
  points( thr.seq, V.seq-2*s_V.seq, type="l", lty=2 )
  points( thr.seq, R.seq          , col="red", type="l" )
  points( thr.seq, R.seq+2*s_R.seq, col="red", type="l", lty=2 )
  points( thr.seq, R.seq-2*s_R.seq, col="red", type="l", lty=2 )
  legend( "bottomright", legend=c("V","R"), col=c("black","red"),
          lty=1, cex=0.75 )
  
  plot( xsub, zsub, xlab="x", ylab="z", asp=1, pch=".", main="" )
  plot  ( thr.seq, pH.sub, type="l",
          main="", xlab="thr", ylab="")
  points( thr.seq, pH.sub+2*s_pH.sub, type="l", lty=2 )
  points( thr.seq, pH.sub-2*s_pH.sub, type="l", lty=2 )
  plot  ( thr.seq, V.sub, col="black", type="l",
          ylim=c(0,max(V.seq,na.rm=T)),
          main="", xlab="thr", ylab="" )
  points( thr.seq, V.sub+2*s_V.sub, type="l", lty=2 )
  points( thr.seq, V.sub-2*s_V.sub, type="l", lty=2 )
  points( thr.seq, R.sub          , col="red", type="l" )
  points( thr.seq, R.sub+2*s_R.sub, col="red", type="l", lty=2 )
  points( thr.seq, R.sub-2*s_R.sub, col="red", type="l", lty=2 )
  legend( "bottomright", legend=c("V","R"), col=c("black","red"),
          lty=1, cex=0.75 )

### Verifying the UQ formulas

  thr <- 1

  PRA.tbl <- t( sapply( 1:length(l_xz.NL), function(d) {
    PRA( l_xz.NL[[d]][,1], l_xz.NL[[d]][,2], thr) } ) )

  PRA.tbl[1:2,]

  slist_pH <- sd( PRA.tbl[,"pH"] )
  slist_V  <- sd( PRA.tbl[,"V" ] )
  slist_R  <- sd( PRA.tbl[,"R" ] )

  par( mfrow=c(2,3) )

  range.pH   <- range( 0, PRA.tbl[,"pH"  ]           )
  range.V    <- range( 0, PRA.tbl[,"V"   ]           )
  range.R    <- range( 0, PRA.tbl[,"R"   ]           )
  range.s_pH <- range( 0, PRA.tbl[,"s_pH"], slist_pH )
  range.s_V  <- range( 0, PRA.tbl[,"s_V" ], slist_V  )
  range.s_R  <- range( 0, PRA.tbl[,"s_R" ], slist_R  )
  
  hist( PRA.tbl[,"pH"], main="p[H]", xlab="", ylab="", xlim=range.pH )
  hist( PRA.tbl[,"V" ], main="V"   , xlab="", ylab="", xlim=range.V  )
  hist( PRA.tbl[,"R" ], main="R"   , xlab="", ylab="", xlim=range.R  )

  hist( PRA.tbl[,"s_pH"], xlim=range.s_pH, main="sd( p[H] )", xlab="", ylab="" )
    abline( v=slist_pH, col="red", lwd=1 )
  hist( PRA.tbl[,"s_V" ], xlim=range.s_V , main="sd( V )"   , xlab="", ylab="" )
    abline( v=slist_V , col="red", lwd=1 )
  hist( PRA.tbl[,"s_R" ], xlim=range.s_R , main="sd( R )"   , xlab="", ylab="" )
    abline( v=slist_R , col="red", lwd=1 )

## Exercises

  xz <- l_xz.NL[[1]] ; x <- xz[,1] ; z <- xz[,2]
  PRA(  x,  z, thr= 1 )
  PRA( -x, -z, thr=-1 )
