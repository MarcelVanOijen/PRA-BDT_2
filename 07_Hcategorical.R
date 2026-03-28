# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 7: Refining the Hazard III: Categorical PRA

  library(fields)
  library(mvtnorm)
  library(plot.matrix)

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
  
  GP.AR <- function( s0, s, x, Sx.s, phi, x0past=0, Sx.t=0, alpha=0 ) {
    ds  <- as.matrix( dist(s) )
    rx  <- exp( -ds/phi )
    ds0 <- sapply(1:length(x),function(i){dist(rbind(s0,s[i,]))})
    r0  <- exp( -ds0/phi )
    m0  <- t(r0) %*% solve(rx) %*% x + alpha * x0past
    V0  <- Sx.s * (1 - t(r0) %*% solve(rx) %*% r0 ) + Sx.t
    return( c( m0=m0, V0=V0 ) ) }

  ns1  <- 8 ; ns2 <- 8 ; nt <- 9
  xst  <- array( NA, dim=c(ns1,ns2,nt) )
  
  Sx.s <- 0.5 ; phi   <- 1
  Sx.t <- 0.5 ; alpha <- 0.5

  set.seed(11)

  for(it in 1:nt) {
    m11 <- if(it==1) { 0 } else { alpha * xst[1,1,it-1] }
    xst[1,1,it] <- rnorm( 1, m11, sqrt(Sx.t) )
    for( s1 in 1:ns1) {
      for( s2 in 1:ns2) {
        if( s1==1 && s2==1) next
        xsti <- xst[,,it]
        s    <- which( !is.na(xsti), arr.ind=T )
        x    <- xsti[s]
        s0   <- c(s1,s2)
        E0   <- if(it==1) { GP.AR( s0, s, x, Sx.s, phi ) }
                else      { GP.AR( s0, s, x, Sx.s, phi,
                                   x0past=xst[s1,s2,it-1], Sx.t, alpha ) }
        xst[s1,s2,it] <- rnorm(1, mean=E0["m0"], sd=sqrt(E0["V0"]) )
      }
    }
  }

  xst <- (xst-min(xst)) / (max(xst)-min(xst))

  ncols <- 5 ; nrows <- ceiling(nt/ncols) ; range.xst <- range(xst,na.rm=T)
  par( mfrow=c(nrows,ncols), mar=c(0,2,0,1), oma=c(0,0,0,0) )
  for(it in 1:nt) {
    if(it < nt) plot( xst[,,it], main="",
                      asp=T, xlab="", ylab="", breaks=range.xst,
                      col=heat.colors(11),
                      axis.col=NULL, axis.row=NULL, key=NULL )
    else        plot( xst[,,it], main="",
                      asp=T, xlab="", ylab="", breaks=range.xst,
                      col=heat.colors(11),
                      axis.col=NULL, axis.row=NULL, fmt.key="%.1f" )
    title( paste0("x(t=",it,")"), line=-2 )
  }

  ncols <- 5 ; nrows <- ceiling( nt / ncols )
  par( mfrow=c(nrows,ncols), mar=c(2,2,1,1), oma=c(0,0,0,0) )

  maxdist <- sqrt( (ns1-1)^2 + (ns2-1)^2 )
  for(it in 1:nt) {
    vgr <- vgram.matrix( xst[,,it], R=maxdist )
    plot( vgr$d, vgr$vgram,
          xlab="", ylab="" ) 
    title( paste0("x(t=",it,")") )
  }

  par( mfrow=c(1,1), mar=c(4,5,0,0) )
  plot( c(xst[,,1:(nt-1)]), c(xst[,,2:nt]), main="",
        xlab="x(t-1)", ylab="x(t)", asp=1, cex.axis=2, cex.lab=2 )
  abline( 0, 1 )

  fz  <- function( x, xpast=1, k=10 ) {
    1 / (1 + exp( -k * (x + xpast/2 - 0.5 ) ) ) }

  zst         <- array( NA, dim=c(ns1,ns2,nt) )
  zst[,,1   ] <- fz( xst[,,1   ] )
  zst[,,2:nt] <- fz( xst[,,2:nt], xst[,,1:(nt-1)] )
  zst         <- (zst-min(zst)) / (max(zst)-min(zst))

  ncols <- 5 ; nrows <- ceiling( nt / ncols )
  par( mfrow=c(nrows,ncols), mar=c(0,2,0,1), oma=c(0,0,0,0) )
  range.zst <- range(zst,na.rm=T) 
  for(it in 1:nt) {
    if(it < nt) plot( zst[,,it], main="",
                      asp=T, xlab="", ylab="", breaks=range.zst,
                      col=rev(terrain.colors(11)),
                      axis.col=NULL, axis.row=NULL, key=NULL )
    else        plot( zst[,,it], main="",
                      asp=T, xlab="", ylab="", breaks=range.zst,
                      col=rev(terrain.colors(11)),
                      axis.col=NULL, axis.row=NULL, fmt.key="%.1f" )
    title( paste0("z(t=",it,")"), line=-2 )
  }

  par( mfrow=c(1,1), mar=c(4,4,0,0), oma=c(0,0,0,0) )
  plot( xst, zst, xlab="x", ylab="z", asp=1, cex.axis=2, cex.lab=2, main="" )

## Single-Category Single-Threshold PRA

  thr.xst <- quantile( xst, pnorm(-1) )

  PRAst <- array( NA, dim=c(ns1,ns2,6) )
  for( s1 in 1:ns1) {
    for( s2 in 1:ns2) { PRAst[s1,s2,] <-
      PRA( x=xst[s1,s2,], z=zst[s1,s2,], thr=thr.xst )
    }
  }

  par( mfrow=c(1,3), mar=c(0,0,0,4), oma=c(0,0,0,0) )

  plot( PRAst[,,1], main="", asp=T, xlab="", ylab="", breaks=c(0,1),
        col=rev(heat.colors(11)), axis.col=NULL , axis.row=NULL, fmt.key="%.1f" )
  title( "p[H]", line=-1 )

  plot( PRAst[,,2], main="", asp=T, xlab="", ylab="", breaks=c(0,1),
        col=rev(heat.colors(11)), axis.col=NULL , axis.row=NULL, fmt.key="%.1f" )
  title( "V", line=-1 )
  
  plot( PRAst[,,3], main="", asp=T, xlab="", ylab="", breaks=c(0,1),
        col=rev(heat.colors(11)), axis.col=NULL , axis.row=NULL, fmt.key="%.1f" )
  title( "R", line=-1 )

## Two-Category Single-Threshold PRA

  PRA.h <- function( x, xpast, z, thr=0 ) {
    n     <- length(z) ; H    <- which(x < thr) ; nH      <- length(H)
    Ez    <- mean( z ) ; Ez_H <- mean( z[H] )   ; Ez_notH <- mean( z[-H] )
    pH    <- nH / n    ; V    <- Ez_notH - Ez_H ; R       <- Ez_notH - Ez
    
    H1    <- which(x < thr & xpast <  thr)      ; nH1     <- length(H1)
    H2    <- which(x < thr & xpast >= thr)      ; nH2     <- length(H2)
    Ez_H1 <- mean( z[H1] )
    Ez_H2 <- mean( z[H2] )
    pH1   <- nH1 / n   ; V1 <- Ez_notH - Ez_H1  ; R1      <- pH1 * V1
    pH2   <- nH2 / n   ; V2 <- Ez_notH - Ez_H2  ; R2      <- pH2 * V2
    return( c(pH=pH,V=V,R=R, pH1=pH1,V1=V1,R1=R1, pH2=pH2,V2=V2,R2=R2) ) }

  PRAst.h <- array( NA, dim=c(ns1,ns2,9) )
  for( s1 in 1:ns1) {
    for( s2 in 1:ns2) { PRAst.h[s1,s2,] <-
      PRA.h( x=xst[s1,s2,2:nt], xpast=xst[s1,s2,1:(nt-1)],
             z=zst[s1,s2,2:nt], thr=thr.xst )
    }
  }

  par( mfrow=c(2,3), mar=c(0,0,0,4), oma=c(0,0,0,0) )

  titles <- c( "p[H]","V","R", "p[H1]","V1","R1", "p[H2]","V2","R2" )
  for( p in 4:8) {
    plot( PRAst.h[,,p], main="", asp=T, xlab="", ylab="", breaks=c(0,1),
          col=rev(heat.colors(11)), axis.col=NULL, axis.row=NULL, key=NULL )
    title( titles[p], line=-1 ) }
  plot(   PRAst.h[,,9], main="", asp=T, xlab="", ylab="", breaks=c(0,1),
          col=rev(heat.colors(11)), axis.col=NULL, axis.row=NULL, fmt.key="%.1f" )
  title( titles[9], line=-1 )
