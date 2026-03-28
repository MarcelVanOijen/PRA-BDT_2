# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 6: Refining the Hazard II: Continuous PRA

  library(mvtnorm)

## Single-Threshold Continuous PRA

  Ez_Gauss <- function( m.=m, S.=S, thr.=thr ) {
    mx     <- m.[1]   ; mz <- m.[2]
    Sx     <- S.[1,1] ; Sz <- S.[2,2]   ; Sxz    <- S.[1,2]
    pthr   <- dnorm(thr., mx, sqrt(Sx)) ; Fthr   <- pnorm(thr., mx, sqrt(Sx))
    Ez_xlo <- mz - Sxz * pthr / Fthr    ; Ez_xhi <- mz + Sxz * pthr / (1-Fthr)
    result <- c( Ez_xlo, Ez_xhi ) ; names ( result ) <- c( "Ez_xlo", "Ez_xhi" )
    return( result ) }

  mx <- mz <- 0          ; m      <- c( mx, mz )
  Sx <- Sz <- 1          ; rxz    <- 0.5
  S  <- diag( c(Sx,Sz) ) ; S[1,2] <- S[2,1] <- rxz * sqrt(Sx * Sz)
  
  px <- function(x, m=mx, S=Sx) {dnorm(x, m, sqrt(S))}
    Ez_x <- function(x){ mz + (x-mx)*rxz }
  v  <- function(x,thr=1) { Ez_Gauss(m,S,thr)["Ez_xhi"] - Ez_x(x) }
  r  <- function(x,thr=1) { px(x) * v(x,thr) }

  x.seq  <- seq( mx-2, mx+2, length.out=41 ) ; thr <- 1
  px.seq <- px( x.seq )
  v.seq  <- v ( x.seq, thr )
  r.seq  <- r ( x.seq, thr )

  par( mfrow=c(1,3), mar=c(3,3,1,0), bty="n", lwd=2 )
  plot( x.seq, px.seq, type="l", xlab="" , xlim=c(-2,1), main="p[x]" )
    abline( v=thr, col="red", lty=2 )
  plot( x.seq, v.seq , type="l", xlab="" , xlim=c(-2,1), main="v(x)" )
    abline( v=thr, col="red", lty=2 )
  plot( x.seq, r.seq , type="l", xlab="x", xlim=c(-2,1), main="r(x)" )
    abline( v=thr, col="red", lty=2 )

## Zero-Threshold Continuous PRA

  par( mfrow=c(1,2), mar=c(4,4,1,4) )

  curve( 1-exp(-x), xlim=c(0,3), ylab="", main="E[z|x]" )

  px_NL <- function(x){ rep(1/3,length(x)) }
  v_NL  <- function(x){ exp(-x) }
  r_NL  <- function(x){ exp(-x) / 3 }

  n_x       <- 16
  x_NL.seq  <- seq( 0, 3, length.out=n_x )
  px_NL.seq <- px_NL( x_NL.seq )
  v_NL.seq  <- v_NL ( x_NL.seq )
  r_NL.seq  <- r_NL ( x_NL.seq )
  v_NLrange <- range( 0, v_NL.seq )
  
  par( mfrow=c(1,3), mar=c(5,3,2,0), bty="n", lwd=2 )
  plot( x_NL.seq, px_NL.seq, type="l", xlab="x", main="p[x]",
        xlim=c(0,3), ylim=c(0,1) )
  plot( x_NL.seq, v_NL.seq , type="l", xlab="x", main="v(x)",
        xlim=c(0,3), ylim=c(0,1) )
  plot( x_NL.seq, r_NL.seq , type="l", xlab="x", main="r(x)",
        xlim=c(0,3), ylim=c(0,1) )

  R.seq  <- ( exp(-x_NL.seq[-n_x]) - exp(-x_NL.seq[-1]) ) / 3
  pH.seq <- ( x_NL.seq[-1] - x_NL.seq[-n_x] ) / 3
  # V.seq  <- R.seq / pH.seq
  V.seq  <- ( exp(-x_NL.seq[-n_x]) - exp(-x_NL.seq[-1])  ) /
            (      x_NL.seq[-1]    -      x_NL.seq[-n_x] )
  
  par( mfrow=c(1,3), mar=c(5,3,2,0) )
  xtxt <- "Upper bound of interval" ; nms <- x_NL.seq[-1]
  barplot( pH.seq, main="p[H]", xlab=xtxt, ylab="", ylim=c(0,0.1), names=nms )
  barplot( V.seq , main="V"   , xlab=xtxt, ylab="", ylim=c(0,1  ), names=nms )
  barplot( R.seq , main="R"   , xlab=xtxt, ylab="", ylim=c(0,0.1), names=nms )

## Loss Distributions
### Gaussian $E[z|x]$

  xz <- rmvnorm( 1e4, m, S ) ; x <- xz[,1] ; z <- xz[,2] ; thr <- 1
  lz <- mean( z[x>=thr] ) - z
  lx <- v(x)
  
  par( mfrow=c(1,2) )

  hist( lz, freq=F, xlim=c(-3,3), ylim=c(0,1), xlab="", main="p[l_z]" )
  lz.seq <- seq(min(lz), max(lz), length = 40)
  lines( lz.seq, dnorm(lz.seq,0.763,1), col="red" )
  
  hist( lx, freq=F, xlim=c(-3,3), ylim=c(0,1), xlab="", main="p[l_x]" )
  lx.seq <- seq(min(lx), max(lx), length = 40)
  lines( lx.seq, dnorm(0.763-lx.seq,0,0.5), col="red" )

### Negative Exponential $E[z|x]$

  x  <- runif( 1e4, 0, 3 )
  z  <- 1 - exp(-x) + rnorm( 1e4, 0, 0.1 ) ; zmax <- 1
  lz <- zmax - z
  lx <- v_NL(x)
  
  par( mfrow=c(1,2) )

  hist( lz, freq=F, xlim=c(-0.5,1.5), ylim=c(0,7), xlab="", main="p[l_z]" )
  lz.seq <- seq(min(lz), max(lz), length = 40)

  hist( lx, freq=F, xlim=c(-0.5,1.5), ylim=c(0,7), xlab="", main="p[l_x]" )
  lx.seq <- seq(min(lx), max(lx), length = 40)
  lines( lx.seq, 1/(3*lx.seq), col="red" )
