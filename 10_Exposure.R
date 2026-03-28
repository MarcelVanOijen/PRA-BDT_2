# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 10: A Third Component of Risk: Exposure

  library(terra)

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

## An Example of Three-Component PRA

  set.seed(1)

  n_t <- 4 ; n_s1 <- 3 ; n_s2 <- 5 ; x <- z <- array( dim=c(n_s1,n_s2,n_t) )
  for(r in 1:n_s1){ for(c in 1:n_s2){ x[r,c,] <- rbeta( n_t, 2, 2 ) } }
  z   <- x^2 + array( rnorm(length(z),sd=0.1), dim=c(n_s1,n_s2,n_t) )

  plot( rast(x), nr=1, reset=T, pax=list(lab=F), main=paste0("x[",1:n_t,"]"),
        range=c(0,1) )

  plot( rast(z), nr=1, reset=T, pax=list(lab=F), main=paste0("z[",1:n_t,"]"),
        range=c(0,1) )

  thr <- 0.25

  whereQ <- function( x=array(dim=c(nlon,nlat,n_t)),
                      z=array(dim=c(nlon,nlat,n_t)), thr.=thr ) {
    ns    <- prod( dim(x)[1:2] ) ; n_t <- dim(x)[3]
    freqH <- function(x,thr.=thr){ sum(x<thr.) }
    n_tH  <- apply(x, c(1,2), freqH)
    siteQ <- which( n_tH > 0, arr.ind=TRUE ) ; nQ <- dim(siteQ)[1]
    return( siteQ ) }

  siteQ <- whereQ( x, z, thr )
  Q     <- matrix( 0, nrow=n_s1, ncol=n_s2 ) ; Q[siteQ] <- 1
  plot( rast(Q), main="Q", grid=T, col=c("green","red"), # R-package 'terra'
        pax=list(xat=1:5,yat=3:1,xlabs=1:5,ylabs=1:3,hadj=2) ) # R-package 'terra'

  PRA3 <- function( x=array(dim=c(nlon,nlat,n_t)),
                    z=array(dim=c(nlon,nlat,n_t)), thr.=thr ) {
    ns    <- prod( dim(x)[1:2] ) ; n_t <- dim(x)[3]
    freqH <- function(x,thr.=thr){ sum(x<thr.) }
    n_tH  <- apply(x, c(1,2), freqH)
    siteQ <- which( n_tH > 0, arr.ind=TRUE ) ; nQ <- dim(siteQ)[1]
    Q     <- nQ / ns ; s_Q <- sqrt( Q*(1-Q) / ns ) 
    
    xQ    <- sapply( 1:n_t, function(i){x[,,i][siteQ]} )
    zQ    <- sapply( 1:n_t, function(i){z[,,i][siteQ]} )
    PRAQ  <- PRA( c(xQ), c(zQ), thr. )
    pH    <- PRAQ["pH"]   ; V   <- PRAQ["V"]   ; R.Q   <- PRAQ["R"]
    s_pH  <- PRAQ["s_pH"] ; s_V <- PRAQ["s_V"] ; s_R.Q <- PRAQ["s_R"]
    R     <- Q * R.Q
    s_R   <- sqrt( s_Q^2*s_R.Q^2 + s_Q^2*R.Q^2 + Q^2*s_R.Q^2 )
    
    result           <- c(  Q ,  pH ,  V ,  R ,  s_Q ,  s_pH ,  s_V ,  s_R  )
    names ( result ) <- c( "Q", "pH", "V", "R", "s_Q", "s_pH", "s_V", "s_R" )
    return( result ) }

  PRA3( x, z, thr )

  PRAst3 <- PRA3( x, z, thr )      ; round(PRAst3,3)

  PRA( c(x), c(z), thr )

  PRAst  <- PRA( c(x), c(z), thr ) ; round(PRAst,3)
