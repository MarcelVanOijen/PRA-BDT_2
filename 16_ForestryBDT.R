# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 16: Example: BDT for Forestry Across Scotland

  library(geodata)
  library(terra)

  set.seed(1)
  
  gadm5  <- world(resolution=5, path="data")
  r_alt  <- rast( "data/r_alt.tif" )
  s_prec <- rast( "data/s_prec.tif" )

  YC <- function( x1, x2 ) {
    alt <- x1  ; prec <- x2
    YC  <- max(0, 10 * (1-exp(-prec/1000)) * (1-alt/1000) )
  return( YC ) }
    
  r_prec_2000 <- s_prec[[1]]
  r_YC_2000   <- xapp( r_alt, r_prec_2000, fun=YC )
  s_YC        <- rast( sapply( 1:19, function(i){
    xapp( r_alt, s_prec[[i]], fun=YC )} ) )

## Multiple Action Levels

  kz <- 30 ; ka <- 0.1
  nI    <- 10 ; layers <- 1:nI ; IRRIG <- (layers-1)*50
  s_u.a <- NULL
  for(i in layers) {
    r_YC.a <- xapp( r_alt, r_prec_2000 + IRRIG[i], fun=YC )
    r_u.a  <- r_YC.a * kz - IRRIG[i] * ka
    s_u.a  <- c( s_u.a, r_u.a ) }
  s_u.a <- rast(s_u.a)
  a.opt <- (which.max( s_u.a ) - 1) * 50

  firrig <- global( a.opt>0 , sum, na.rm=T ) / global( a.opt>=0, sum, na.rm=T )

  par( mfrow=c(1,2), mar=c(2,2,2,3), cex=0.8 )

  plot.aopt <- plot( a.opt, main=expression( 'Optimum irrigation (mm y'^-1*')' ),
        xaxt="n", yaxt="n",
        col=terrain.colors(100,rev=T)) ; plot( gadm5, add=T)

  global( a.opt>0, sum, na.rm=T ) / global( a.opt>=0, sum, na.rm=T )

## PRA vs. BDT

  IRRIG      <- 500
  s_YC       <- rast( sapply(1:19, function(i){
    xapp( r_alt, s_prec[[i]]      , fun=YC )}) )
  s_YC.IRRIG <- rast( sapply(1:19, function(i){
    xapp( r_alt, s_prec[[i]]+IRRIG, fun=YC )}) )

  kz        <- 30 ; ka <- 0.1
  s_u       <- s_YC       * kz
  s_u.IRRIG <- s_YC.IRRIG * kz - IRRIG * ka

  par( mfrow=c(2,4), mar=c(1,0,1,0) )

  years <- c(1:2,NA,19)

  for( y in years ) { 
    if(is.na(y)) plot(-1:1,rep(0,3),axes=0,xlim=c(-5,5),pch=19)
    else { plot( s_u[[y]], range=c(0,230),
                 main=paste0(y+1999,"\nIRRIG=0"), xaxt='n', yaxt='n' )
           plot( gadm5, add=T) } }

  for( y in years ) { 
    if(is.na(y)) plot(-1:1,rep(0,3),axes=0,xlim=c(-5,5),pch=19)
    else { plot( s_u.IRRIG[[y]], range=c(0,230),
                 main=paste0(y+1999,"\nIRRIG=",IRRIG), xaxt='n', yaxt='n' )
           plot( gadm5, add=T) } }

  r_u       <- mean( s_u )
  r_u.IRRIG <- mean( s_u.IRRIG )

  par( mfrow=c(1,3), mar=c(1,1,1,3) )

  plot( r_u            , main="E[ u | IRRIG=0 ]",
        xaxt='n', yaxt='n', range=c(  0,230) ) ; plot( gadm5, add=T)
  plot( r_u.IRRIG      , main=paste0("E[ u | IRRIG=",IRRIG," ]"),
        xaxt='n', yaxt='n', range=c(  0,230) ) ; plot( gadm5, add=T)
  plot( r_u.IRRIG - r_u, main="d(u)",
        xaxt='n', yaxt='n', range=c(-50, 10) ) ; plot( gadm5, add=T)

  PRA.modif <- function( x, z, thr=1000, rel=F ) {
    n       <- length(x) ; H  <- which(x < thr) ; notH   <- which(x >= thr)
    n_H     <- length(H) ; pH <- n_H / n        ; n_notH <- length(notH)
    Ez_H    <- mean( z[H]    ) ; s_Ez_H    <- sqrt( var(z[H   ]) / n_H    )
    Ez_notH <- mean( z[notH] ) ; s_Ez_notH <- sqrt( var(z[notH]) / n_notH )
    if(n_H==0){ Ez_H <- min(z) } ; if(n_notH==0){ Ez_notH <- max(z) }
    V       <- Ez_notH - Ez_H  ; R         <- pH * V
    s_pH    <- sqrt( pH*(1-pH) / n )
    s_V     <- sqrt( s_Ez_H^2 + s_Ez_notH^2 )
    s_R     <- sqrt( s_pH^2*s_V^2 + s_pH^2*V^2 + pH^2*s_V^2 )
    if(rel){V <- V/Ez_notH ; R=R/Ez_notH; s_V=s_V/Ez_notH; s_R=s_R/Ez_notH}
    return( c(Ez_notH=Ez_notH, pH=pH, V=V, R=R, s_pH=s_pH, s_V=s_V, s_R=s_R) )
  }

  s_PRA       <- xapp( s_prec      , s_u      , fun=PRA.modif )
  s_PRA.IRRIG <- xapp( s_prec+IRRIG, s_u.IRRIG, fun=PRA.modif )
  r_Ez_notH   <- s_PRA$Ez_notH    ; r_Ez_notH.IRRIG <- s_PRA.IRRIG$Ez_notH
  r_pH        <- s_PRA$pH         ; r_pH.IRRIG      <- s_PRA.IRRIG$pH
  r_pH        <- mask(r_pH,r_alt) ; r_pH.IRRIG      <- mask(r_pH.IRRIG,r_alt)
  r_V         <- s_PRA$V          ; r_V.IRRIG       <- s_PRA.IRRIG$V
  r_R         <- s_PRA$R          ; r_R.IRRIG       <- s_PRA.IRRIG$R

  par( mfrow=c(2,3), mar=c(1,1,1,3) )

  plot( r_pH      , main="p[H]", xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0, 1) ) ; plot( gadm5, add=T)
  plot( r_V       , main="V" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,32) ) ; plot( gadm5, add=T)
  plot( r_R       , main="R" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,23) ) ; plot( gadm5, add=T)

  plot( r_pH.IRRIG, main="p[H].IRRIG", xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0, 1) ) ; plot( gadm5, add=T)
  plot( r_V.IRRIG , main="V.IRRIG" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,32) ) ; plot( gadm5, add=T)
  plot( r_R.IRRIG , main="R.IRRIG" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T), range=c( 0,23) ) ; plot( gadm5, add=T)

  r_Rc       <- r_R       - r_Ez_notH
  r_Rc.IRRIG <- r_R.IRRIG - r_Ez_notH.IRRIG

  par( mfrow=c(1,3), mar=c(1,1,1,3) )

  plot( r_Rc             , main=expression( 'R'[c] ),
        xaxt='n', yaxt='n', range=c(-270, 0) ) ; plot( gadm5, add=T)
  plot( r_Rc.IRRIG       , main=expression( 'R'[c.IRRIG] ),
        xaxt='n', yaxt='n', range=c(-270, 0) ) ; plot( gadm5, add=T)
  plot( r_Rc - r_Rc.IRRIG, main=expression( '-d(R'[c]*')' ),
        xaxt='n', yaxt='n', range=c( -50,10) ) ; plot( gadm5, add=T)
