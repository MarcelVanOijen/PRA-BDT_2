# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 8: Refining the Hazard IV: Multi-Hazard PRA

  library(fields)
  library(mvtnorm)
  library(plot.matrix)

  PRAi <- function( xz, thr=c(0,0) ) {
    x1  <- xz[,1]  ; x2 <- xz[,2] ; z <- xz[,3]
    n_c <- 2^2 - 1 ; n  <- length(x1)
    H   <- vector("list",n_c)
    n_H <- pH <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_c)
    H[[1]]  <- which(x1 <  thr[1] & x2 <  thr[2]) ; n_H[1] <- length(H[[1]])
    H[[2]]  <- which(x1 <  thr[1] & x2 >= thr[2]) ; n_H[2] <- length(H[[2]])
    H[[3]]  <- which(x1 >= thr[1] & x2 <  thr[2]) ; n_H[3] <- length(H[[3]])
    NotH    <- which(x1 >= thr[1] & x2 >= thr[2]) ; n_NotH <- length(NotH)
    pH      <- n_H / n         ; s_pH      <- sqrt( pH*(1-pH) / n )
    Ez_NotH <- mean( z[NotH] ) ; s_Ez_NotH <- sqrt( var(z[NotH] ) / n_NotH )   
    for(i in 1:n_c) {
      Ez_Hi   <-       mean( z[ H[[i]] ] )
      s_Ez_Hi <- sqrt( var ( z[ H[[i]] ] ) / n_H[i] )
      V[i]    <- Ez_NotH - Ez_Hi
      s_V[i]  <- sqrt( s_Ez_NotH^2 + s_Ez_Hi^2 ) }
    R     <- pH * V
    s_R   <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2  )
    return( cbind( 1:3, pH, V, R, s_pH, s_V, s_R ) ) }

  set.seed(1)

  m_G3  <- rep(0,3) ; S_G3  <- matrix(0.5,nrow=3,ncol=3) ; diag(S_G3) <- 1
  n_G3  <- 300      ; xz_G3 <- rmvnorm( n_G3, m_G3, S_G3 )

  panel.hist <- function(x, ...) {
    h      <- hist(x, plot=F)
    breaks <- h$breaks; n_b <- length(breaks)
    y      <- h$counts; y   <- min(x) + (max(x)-min(x)) * y/max(y)
    rect(breaks[-n_b], min(x), breaks[-1], y, col="cyan", ...) }

  pairs(xz_G3, diag.panel=panel.hist, label="" )

  PRAi( xz_G3, c(0,0) )

## Dataset with Two Uncorrelated Hazard Variables that Exert Correlated V

  set.seed(1)

  n  <- 3e2
  x1 <- rbeta( n, 3, 3 ) ; x2 <- rbeta( n, 3, 3 )
  z  <- as.integer( x1 >= 0.5 | x2 >= 0.5 )
  xz <- cbind( x1, x2, z )

  pairs( xz, diag.panel=panel.hist, label="" )

  PRAi( xz, thr=c(0.5,0.5) )
