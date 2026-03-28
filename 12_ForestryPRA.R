# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 12: Example: PRA for Forestry Across Scotland

  library(geodata)
  library(terra)

  set.seed(1)
  
## Data
  
  spdf_GBR <- gadm( country="GBR", level=0, path="data" ) # Great Britain
  spdf_SCO <- crop( spdf_GBR, ext(-7.7,-1.8,54.5,58.6) ) # Scotland
  ext_SCO  <- ext( spdf_SCO )
  gadm5    <- world(resolution=5, path="data")
  
  # dir_WC_prec.2000_2009 <- paste0( dir_WC, "prec_2000_2009/" )
  # dir_WC_prec.2010_2018 <- paste0( dir_WC, "prec_2010_2018/" )
  # s_prec.mon            <- rast( c(
  #   list.files( dir_WC_prec.2000_2009, ".tif", full.names=T ),
  #   list.files( dir_WC_prec.2010_2018, ".tif", full.names=T ) ) )
  # s_prec.mon            <- crop(s_prec.mon, ext_SCO)
  # s_prec.mon            <- s_prec.mon[[1:228]] # Jan 2000 to Dec 2018
  # names(s_prec.mon)     <- stringr::str_sub( names(s_prec.mon), -7, -1 )
  # years                 <- rep( 2000:2018, each=12 )
  # s_prec                <- tapp( s_prec.mon, index=years, fun=sum )
  # writeRaster( s_prec, "data/climate/s_prec.tif", overwrite=T )
  s_prec      <- rast( "data/s_prec.tif" )
  
  # r_alt_GBR   <- elevation_30s( country='GBR', path='data', mask=F )
  # r_alt.hires <- crop( r_alt_GBR, ext_SCO )
  # r_alt       <- resample( r_alt.hires, s_prec )
  # r_tree      <- landcover(var='trees', path='data', download=F)
  # r_tree      <- crop     ( r_tree, ext_SCO )
  # r_tree      <- mask     ( r_tree, r_alt.hires )
  # writeRaster( r_alt , "data/elevation/r_alt.tif", overwrite=T )
  # writeRaster( r_tree, "data/landuse/r_tree.tif" , overwrite=T )
  r_alt  <- rast( "data/r_alt.tif" )
  r_tree <- rast( "data/r_tree.tif" )
  r_tree <- resample ( r_tree, s_prec )

  par( mfrow=c(1,3), mar=c(2,2,2,3) )

  plot( r_alt, col=terrain.colors(100), main="Altitude (m)" )
        plot( gadm5, add=T )
  plot( s_prec[[ 1]], main="Prec. in 2000 (mm)",
        col=rainbow(100,end=0.7), range=c(0,3000), xaxt="n", yaxt="n" )
        plot( gadm5, add=T )
  plot( s_prec[[19]], main="Prec. in 2018 (mm)",
        col=rainbow(100,end=0.7), range=c(0,3000), xaxt="n", yaxt="n" )
        plot( gadm5, add=T )

## Two-component Model-Based PRA

  YC <- function(alt,prec){max(0, 10 * (1-exp(-prec/1000)) * (1-alt/1000))}

  r_prec_2000 <- s_prec[[1]]
  r_YC_2000   <- xapp( r_alt, r_prec_2000, fun=YC )
  s_YC        <- rast( sapply( 1:19, function(i){
    xapp( r_alt, s_prec[[i]], fun=YC )} ) )

  par( mfrow=c(1,2), mar=c(2,2,2,3) )

  plot( r_YC_2000, main="Yield class", col=terrain.colors(100,rev=T),
        range=c(0,9), plg=list(cex=.7), pax=list(side=1:2, cex=.7) )
        plot( gadm5, add=T )
  plot( r_tree, main="Tree cover", col=terrain.colors(100,rev=T),
        plg=list(cex=.7), pax=list(side=1:2, cex=.7) )
        plot( gadm5, add=T )

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

  s_PRA <- xapp( s_prec, s_YC, fun=PRA.modif )
  r_pH  <- s_PRA$pH ; r_V <- s_PRA$V ; r_R <- s_PRA$R

  r_pH  <- mask(r_pH,r_alt)

  par( mfrow=c(1,3), mar=c(2,2,2,3) )

  plot( r_pH, main="p[H]", xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )
  plot( r_V , main="V" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )
  plot( r_R , main="R" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )

  r_s_pH  <- s_PRA$s_pH         ; r_s_V <- s_PRA$s_V ; r_s_R <- s_PRA$s_R
  r_s_pH  <- mask(r_s_pH,r_alt)

  par( mfrow=c(1,3), mar=c(2,2,2,3) )

  plot( r_s_pH, main="s_pH", xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )
  plot( r_s_V , main="s_V" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )
  plot( r_s_R , main="s_R" , xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )

  cor( terra::values(r_pH ), terra::values(r_R), use="pairwise.complete.obs" )^2
  cor( terra::values(r_V  ), terra::values(r_R), use="pairwise.complete.obs" )^2

## Three-Component Model-Based PRA

  r_Q  <- r_tree
  r_Rg <- r_R * r_Q

  par( mfrow=c(1,2), mar=c(3,1,1,5) )

  plot( r_Rg, main=expression( 'R'[g] ), xaxt='n',yaxt='n',
        col=terrain.colors(100,rev=T) ) ; plot( gadm5, add=T )
  hist( r_Rg, main=expression( 'R'[g] ) )

  cor( terra::values(r_pH ), terra::values(r_Rg), use="pairwise.complete.obs" )^2
  cor( terra::values(r_Q  ), terra::values(r_Rg), use="pairwise.complete.obs" )^2
  cor( terra::values(r_V  ), terra::values(r_Rg), use="pairwise.complete.obs" )^2
