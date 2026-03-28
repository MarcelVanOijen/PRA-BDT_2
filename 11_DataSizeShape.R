# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 11: The Impact of Dataset Size and Shape

  library(MCMCpack)
  library(mvtnorm)

  PRA0 <- function( x, z, thr=0 ) {
    n    <- length(x)      ;  H       <- which(x < thr)  ;  n_H <-length(H)
    Ez_H <- mean( z[ H] )  ;  Ez_notH <- mean( z[-H] )
    pH   <- n_H / n        ;  V       <- Ez_notH - Ez_H  ;  R   <- pH * V
    return( c( pH=pH, V=V, R=R ) ) }

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

  PRA0_Gauss <- function( m.=m, S.=S, thr.=thr ) {
    mx <- m.[1] ; sx <- sqrt(S.[1,1]) ; Sxz <- S.[1,2]
    pH <- pnorm( thr., mx, sx )
    V  <- Sxz * dnorm(thr., mx, sx) / (pH * (1-pH))
    R  <- pH * V
    return( c( pH=pH, V=V, R=R ) ) }

  PRA_Gauss <- function( m.=m, S.=S, n.=n, thr.=thr ) {
    pH <- V <- R <- rep( NA, 1e3 )
    for(j in 1:1e3){
      S     <- riwish( n.-1, S. * (n.-1) ) ; m <- rmvnorm( 1, m., S/n. )
      PRA   <- PRA0_Gauss( m, S, thr. )
      pH[j] <- PRA["pH"] ; V[j] <- PRA["V"] ; R[j] <- PRA["R"]
    }
    return( c( pH=mean(pH),   V=mean(V),   R=mean(R),
             s_pH=sd  (pH), s_V=sd  (V), s_R=sd  (R) ) )
  }

  PRA_BayesL <- function( x, z, thr=0, se=1 ) {
    n       <- length(x)  ;  X  <- cbind(1,x)
    mb      <- c(0,0)     ;  Sb <- diag(1.e4,2)  ;  Se <- diag(se^2,n)
    Sb_y    <- solve( solve(Sb) + t(X) %*% solve(Se) %*% X )
    mb_y    <- Sb_y %*% (solve(Sb) %*% mb + t(X) %*% solve(Se) %*% z)
    i_H     <- which( x <  thr )       ; n_H     <- length(i_H)
    i_NotH  <- which( x >= thr )       ; n_NotH  <- length(i_NotH)
    pH      <- n_H / n                 ; s_pH    <- sqrt( pH*(1-pH) / n )
    Ex_H    <- mean( x[i_H] )          ; Ex_NotH <- mean( x[i_NotH] )
    Ez_H    <- c(1,Ex_H) %*% mb_y ; Ez_NotH <- c(1,Ex_NotH) %*% mb_y
    V       <- Ez_NotH - Ez_H          ; R       <- pH * V
      Szi   <- function(i){ t(c(1,x[i])) %*% Sb_y %*% c(1,x[i]) + se^2 }
    Sz_H    <- sum( sapply(i_H   ,Szi) ) / n_H    + mb_y[2]^2 * var(x[i_H])
    Sz_NotH <- sum( sapply(i_NotH,Szi) ) / n_NotH + mb_y[2]^2 * var(x[i_NotH])
    s_V     <- sqrt( Sz_H / n_H + Sz_NotH / n_NotH )
    s_R     <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2 )
    return( list( mb  = mb_y, Sb = Sb_y,
                  PRA = c( pH=pH, V=V, R=R, s_pH=s_pH, s_V=s_V, s_R=s_R ) ) )
  }

## Sparse Dataset ($n=4$)

  load( "data/xz_4.RData" )
  m_4 <- c( mean(x_4), mean(z_4) ) ; S_4 <- cov( cbind(x_4,z_4) ) ; n_4 <- 4

  par( mfrow=c(1,1) )
  plot( x_4, z_4, xlab="x", ylab="z", bty="n", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5) )

  set.seed(1)

  PRA_Sampling     <- PRA       ( x_4, z_4,      0 )
  PRA_Distribution <- PRA_Gauss ( m_4, S_4, n_4, 0 )
  PRA_Model        <- PRA_BayesL( x_4, x_4,      0 )$PRA

  rbind( PRA_Sampling, PRA_Distribution, PRA_Model )

## A Sequence of Linear Datasets ($n$ Exponentially Increasing)

  set.seed(1)

  mu  <- c(0,0) ; Sigma  <- diag(1,2) ; Sigma[1,2] <- Sigma[2,1] <- 0.5
  n_d <- 9      ; l_xz.Ls <- vector("list",n_d)
  for(d in 1:n_d) { l_xz.Ls[[d]] <- rmvnorm( 2^(d+2), mu, Sigma ) }

  par( mfrow=c(1,3), mar=c(5,4,2,4) )
  plot( l_xz.Ls[[1  ]][,1], l_xz.Ls[[1  ]][,2], main="dataset 1",
        xlab="x", ylab="z" )
  plot( l_xz.Ls[[5  ]][,1], l_xz.Ls[[5  ]][,2], main="dataset 5",
        xlab="x", ylab="z" )
  plot( l_xz.Ls[[n_d]][,1], l_xz.Ls[[n_d]][,2], main=paste("dataset",n_d),
        xlab="x", ylab="z" )

  thr     <- 0
  PRA.tbl <- t( sapply( 1:n_d, function(d){
    PRA ( l_xz.Ls[[d]][,1], l_xz.Ls[[d]][,2], thr=thr) } ) )

  PRA.tbl2 <- t( sapply( 1:n_d, function(d){
    PRA_Gauss( m.=colMeans(l_xz.Ls[[d]]), cov(l_xz.Ls[[d]]), 2^(d+2), thr) } ) )
  PRA.tbl3 <- t( sapply( 1:n_d, function(d){
    PRA_BayesL( l_xz.Ls[[d]][,1], l_xz.Ls[[d]][,2], thr=0, se=1 )$PRA } ) )

  r <- 1 ; n <- rep(2^(r+2),3)
  results.lo  <- cbind( n, rbind( PRA.tbl[r,], PRA.tbl2[r,], PRA.tbl3[r,] ) )
  rownames(results.lo) <- c("Sampling-based", "Distribution-based", "Model-based" )
  round( results.lo, 3 )

  r <- 9 ; n <- rep(2^(r+2),3)
  results.hi <- cbind( n, rbind( PRA.tbl[r,], PRA.tbl2[r,], PRA.tbl3[r,] ) )
  rownames(results.hi) <- c("Sampling-based", "Distribution-based", "Model-based" )
  round( results.hi, 3 )

  par(mfrow=c(2,3), mar=c(5,2,2,1))

  plot  ( 1:n_d, PRA.tbl [,"pH"], main="pH", xlab="", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(3:11) )
  points( 1:n_d, PRA.tbl2[,"pH"], pch=20, col="blue" )
  points( 1:n_d, PRA.tbl3[,"pH"], pch= 3, col="green" )
  abline(h=0.5  ,col="red",lty=2)
    
  plot  ( 1:n_d, PRA.tbl [,"V" ], main="V" , xlab="", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(3:11) )
  points( 1:n_d, PRA.tbl2[,"V" ], pch=20, col="blue" )
  points( 1:n_d, PRA.tbl3[,"V" ], pch= 3, col="green" )
  abline(h=0.798,col="red",lty=2)

  plot  ( 1:n_d, PRA.tbl [,"R" ], main="R" , xlab="", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(3:11) )
  points( 1:n_d, PRA.tbl2[,"R" ], pch=20, col="blue" )
  points( 1:n_d, PRA.tbl3[,"R" ], pch= 3, col="green" )
  abline(h=0.399,col="red",lty=2)
  
  legend ( "bottomright", col=c("black","blue","green"), pch=c(1,20,3), cex=0.75,
           legend=c("Sampling-based","Distribution-based","Model-based") )

  plot  ( 1:n_d, PRA.tbl [,"s_pH"], main="s_pH", xlab="n", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(3:11) )
  points( 1:n_d, PRA.tbl2[,"s_pH"], pch=20, col="blue" )
  points( 1:n_d, PRA.tbl3[,"s_pH"], pch= 3, col="green" )

  plot  ( 1:n_d, PRA.tbl [,"s_V" ] , main="s_V" , xlab="n", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(3:11) )
  points( 1:n_d, PRA.tbl2[,"s_V" ], pch=20, col="blue" )
  points( 1:n_d, PRA.tbl3[,"s_V" ], pch= 3, col="green" )

  plot  ( 1:n_d, PRA.tbl [,"s_R" ] , main="s_R" , xlab="n", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(3:11) )
  points( 1:n_d, PRA.tbl2[,"s_R" ], pch=20, col="blue" )
  points( 1:n_d, PRA.tbl3[,"s_R" ], pch= 3, col="green" )

  legend ( "topright", col=c("black","blue","green"), pch=c(1,20,3), cex=0.75,
           legend=c("Sampling-based","Distribution-based","Model-based") )

## A Sequence of Nonlinear Datasets ($n$ Exponentially Increasing)

  set.seed(1)

  n_d <- 10 ; l_xz.NLs <- vector("list",n_d) ; sz <- 0.1
  for(d in 1:n_d) {
    x   <- runif( 2^(d+1), 0, 3 )
    ez  <- rnorm( 2^(d+1), 0, sz) ; z <- 1-exp(-x) + ez
    l_xz.NLs[[d]] <- cbind(x,z) }

  par( mfrow=c(1,3), mar=c(5,4,2,4) )
  plot( l_xz.NLs[[1  ]][,1], l_xz.NLs[[1  ]][,2], main="dataset 1",
        xlab="x", ylab="z", ylim=c(-0.5,1.5) )
  plot( l_xz.NLs[[6  ]][,1], l_xz.NLs[[6  ]][,2], main="dataset 6",
        xlab="x", ylab="z", ylim=c(-0.5,1.5) )
  plot( l_xz.NLs[[n_d]][,1], l_xz.NLs[[n_d]][,2], main=paste("dataset",n_d),
        xlab="x", ylab="z", ylim=c(-0.5,1.5) )

  thr <- 1
  PRA.tbl <- t( sapply( 1:n_d, function(d){
    PRA( l_xz.NLs[[d]][,1], l_xz.NLs[[d]][,2], thr=thr ) } ) )
  PRA.tbl2 <- t( sapply( 1:n_d, function(d){
    PRA_Gauss( m.=colMeans(l_xz.NLs[[d]]), cov(l_xz.NLs[[d]]), 2^(d+1), thr ) } ) )
  # PRA.tbl3 <- t( sapply( 1:n_d, function(d){
  PRA.tbl3 <- t( sapply( c(1,3:n_d), function(d){ # No values of x<1 in dataset 2
    PRA_BayesL( l_xz.NLs[[d]][,1], l_xz.NLs[[d]][,2], thr=thr, se=1 )$PRA } ) )

  n <- rep(4,3)
  results.lo <- cbind( n, rbind( PRA.tbl[1,], PRA.tbl2[1,], PRA.tbl3[1,] ) )
  rownames(results.lo) <- c("Sampling-based", "Distribution-based", "Model-based" )
  round( results.lo, 3 )

  n <- rep(2048,3)
  results.hi <- cbind( n, rbind( PRA.tbl[10,], PRA.tbl2[10,], PRA.tbl3[9,] ) )
  rownames(results.hi) <- c("Sampling-based", "Distribution-based", "Model-based" )
  round( results.hi, 3 )

  par(mfrow=c(2,3), mar=c(5,2,2,1))

  PRA_inf <- PRA( l_xz.NLs[[d]][,1], l_xz.NLs[[d]][,2], thr=thr )

  ylim_pH <- range( 0, PRA.tbl[,"pH"], PRA.tbl2[,"pH"], PRA.tbl3[,"pH"] )
  ylim_V  <- range( 0, PRA.tbl[,"V" ], PRA.tbl2[,"V" ], PRA.tbl3[,"V" ], na.rm=T )
  ylim_R  <- range( 0, PRA.tbl[,"R" ], PRA.tbl2[,"R" ], PRA.tbl3[,"R" ], na.rm=T )

  seq3 <- c(1,3:n_d)
  
  plot  ( 1:n_d, PRA.tbl [,"pH"], main="pH", xlab="", ylab="", ylim=ylim_pH,
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(2:11) )
  points( 1:n_d, PRA.tbl2[,"pH"], pch=20, col="blue"  )
  points( seq3 , PRA.tbl3[,"pH"], pch=20, col="green" )
    abline(h=PRA_inf["pH"],col="red",lty=2)
    
  plot  ( 1:n_d, PRA.tbl [,"V" ], main="V" , xlab="", ylab="", ylim=ylim_V,
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(2:11) )
  points( 1:n_d, PRA.tbl2[,"V" ], pch=20, col="blue"  )
  points( seq3 , PRA.tbl3[,"V" ], pch=20, col="green" )
    abline(h=PRA_inf["V"] ,col="red",lty=2)
    
  plot  ( 1:n_d, PRA.tbl [,"R" ], main="R" , xlab="", ylab="", ylim=ylim_R,
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(2:11) )
  points( 1:n_d, PRA.tbl2[,"R" ], pch=20, col="blue"  )
  points( seq3 , PRA.tbl3[,"R" ], pch=20, col="green" )
    abline(h=PRA_inf["R"] ,col="red",lty=2)
    legend ( "bottomright", col=c("black","blue","green"), pch=c(1,20), cex=0.75,
             legend=c("Sampling-based","Distribution-based","Model-based") )
    
  ylim_s_pH <- range( 0, PRA.tbl[,"s_pH"], PRA.tbl2[,"s_pH"], PRA.tbl3[,"s_pH"] )
  ylim_s_V  <- range( 0, PRA.tbl[,"s_V" ], PRA.tbl2[,"s_V" ], PRA.tbl3[,"s_V" ], na.rm=T )
  ylim_s_R  <- range( 0, PRA.tbl[,"s_R" ], PRA.tbl2[,"s_R" ], PRA.tbl3[,"s_R" ], na.rm=T )

  plot  ( 1:n_d, PRA.tbl [,"s_pH"], main="s_pH", xlab="n", ylab="",
          xaxt="n" )
  axis  ( 1, at=1:n_d, labels=2^(2:11) )
  points( 1:n_d, PRA.tbl2[,"s_pH"], pch=20, col="blue"  )
  points( seq3 , PRA.tbl3[,"s_pH"], pch=20, col="green" )
  
  plot  ( 1:n_d, PRA.tbl [,"s_V" ], main="s_V" , xlab="n", ylab="",
          xaxt="n", ylim=ylim_s_V )
  axis  ( 1, at=1:n_d, labels=2^(2:11) )
  points( 1:n_d, PRA.tbl2[,"s_V" ], pch=20, col="blue" )
  points( seq3 , PRA.tbl3[,"s_V" ], pch=20, col="green" )
  
  plot  ( 1:n_d, PRA.tbl[,"s_R"  ], main="s_R" , xlab="n", ylab="" ,
          xaxt="n", ylim=ylim_s_V )
  axis  ( 1, at=1:n_d, labels=2^(2:11) )
  points( 1:n_d, PRA.tbl2[,"s_R" ], pch=20, col="blue" )
  points( seq3 , PRA.tbl3[,"s_R" ], pch=20, col="green" )
  legend ( "topright", col=c("black","blue","green"), pch=c(1,20), cex=0.75,
           legend=c("Sampling-based","Distribution-based","Model-based") )
