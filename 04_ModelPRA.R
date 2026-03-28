# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 4: Model-Based PRA

  library(mvtnorm)
  library(nimble)
  library(scales)

  load( "data/l_xz.L.RData" )
  load( "data/l_xz.NL.RData" )

## Example I: Linear Relationship

  xz <- l_xz.L[[1]] ; m  <- colMeans(xz)  
  x  <- xz[,1]      ; X  <- cbind( 1, x ) ; n  <- length(x)
  z  <- xz[,2]      ; se <- 1             ; Se <- diag(se^2,n)

### Model and Conjugate Bayesian Calibration
  
  mb   <- c(0,0)   ;   Sb <- diag(1.e4,2)
  Sb_y <- solve( solve(Sb) + t(X) %*% solve(Se) %*% X )
  mb_y <- Sb_y %*% (solve(Sb) %*% mb + t(X) %*% solve(Se) %*% z)

  smpl_b <- rmvnorm( 10, mb_y, Sb_y )

  par(mfrow=c(1,2))
  plot( x, z ) ; abline( mb_y, col="red" )
  nsmpl <- 10  ; smpl_b <- rmvnorm( nsmpl, mean=mb_y, sigma=Sb_y )
  for( i in 1:nsmpl) { abline( smpl_b[i,], lty=2 ) }

### PRA

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

  PRA_BayesL( x, z )$PRA

  thr <- seq( -2, 2, by=0.05 ) ; n_thr <- length(thr)
  n_H <- n_notH <- Ez_H <- Ez_notH <- Sz_H <- Sz_notH <- numeric(n_thr)
  R   <- V      <- pH   <- s_R     <- s_V  <- s_pH    <- numeric(n_thr)
  LCIR     <- UCIR     <- LCIV  <- UCIV               <- numeric(n_thr)
  LCIz_notH <- UCIz_notH <- LCIz_H <- UCIz_H              <- numeric(n_thr)
  for(i in 1:n_thr) {
    i_H    <- which( x <  thr[i] ) ; n_H[i]    <- length(i_H)
    i_notH <- which( x >= thr[i] ) ; n_notH[i] <- length(i_notH)
    pH[i]  <- n_H[i] / n           ; s_pH[i]   <- sqrt( pH[i]*(1-pH[i]) / n )
    Ex_H       <- mean( x[i_H] )    ; Ez_H[i]    <- c(1,Ex_H)    %*% mb_y
    Ex_notH    <- mean( x[i_notH] ) ; Ez_notH[i] <- c(1,Ex_notH) %*% mb_y
    V[i]       <- Ez_notH[i] - Ez_H[i]
    R[i]       <- pH[i] * V[i]
    Szi        <- function(i){ t(c(1,x[i])) %*% Sb_y %*% c(1,x[i]) + se^2 }
    Sz_H[i]    <- sum( sapply(i_H   ,Szi) ) / n_H[i]    + mb_y[2]^2 * var(x[i_H])
    Sz_notH[i] <- sum( sapply(i_notH,Szi) ) / n_notH[i] + mb_y[2]^2 * var(x[i_notH])
    s_V[i]     <- sqrt( Sz_H[i] / n_H[i] + Sz_notH[i] / n_notH[i] )
    s_R[i]     <- sqrt( s_pH[i]^2 * s_V[i]^2 + s_pH[i]^2 * V[i]^2 + pH[i]^2 * s_V[i]^2 )
  }
  LCIz_notH <- Ez_notH - sqrt( Sz_notH / n_notH )
  UCIz_notH <- Ez_notH + sqrt( Sz_notH / n_notH )
  LCIz_H    <- Ez_H    - sqrt( Sz_H    / n_H    )
  UCIz_H    <- Ez_H    + sqrt( Sz_H    / n_H    )
  LCIV      <- V - s_V ; UCIV <- V + s_V ; LCIR <- R - s_R ; UCIR <- R + s_R 

  par    ( mfrow=c(1,2), mar=c(5,5,1,1) )

  plot   ( thr, Ez_H, type="l",
           ylim=range(c(LCIz_notH,LCIz_H,UCIz_notH,UCIz_H),na.rm=T),
           xlab="Threshold", ylab="Expectation of z")
  polygon( c(thr,thr[n_thr:1]), c(LCIz_H,UCIz_H[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIz_notH,UCIz_notH[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  abline ( h=m[2], col="blue", lty=3 )
  lines  ( thr, Ez_H, col="red") ; lines( thr, Ez_notH, col="black" )
  legend ( "bottomright", legend=c("E[z|¬H]","E[z]","E[z|H]"),
           col=c("black","blue","red"), lty=c(1,3,1), cex=0.75 )

  plot   ( thr, V, col="black", type="l",
           ylim=range(c(LCIV,UCIV,LCIR,UCIR),na.rm=T),
           xlab="Threshold", ylab="V, R" )
  polygon( c(thr,thr[n_thr:1]), c(LCIV,UCIV[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIR,UCIR[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  lines  ( thr, V, col="black" ) ; lines( thr, R, col="red" )
  legend ( "bottomright", legend=c("V","R"), col=c("black","red"), lty=1, cex=0.75 )

## Example II: Nonlinear Relationship

  xz <- l_xz.NL[[1]] ; m <- colMeans(xz) ; S <- cov(xz)
  x  <- xz[,1]       ; z <- xz[,2]       ; n <- length(x)

### Model and MCMC

  Model1.Code <- nimbleCode( {
    alpha  ~ dnorm( 0, sd=100 )
    beta   ~ dnorm( 0, sd=100 )
    tau    ~ dgamma( 0.01, 0.01 )
    sigma <- 1 / sqrt(tau)
    for(i in 1:ndata){
      mu[i] <- alpha + beta*exp(-x[i])
      z[i]   ~ dnorm( mu[i], sd=sigma ) }
  } )

  Model1.Constants <- list( ndata=n, x=x )
  Model1.Data      <- list(z=z)
  Model1.Nimble    <- nimbleModel  ( Model1.Code, constants=Model1.Constants,
                                                  data=Model1.Data )
  Model1.Comp      <- compileNimble( Model1.Nimble )
  Model1.Conf      <- configureMCMC( Model1.Nimble, print=F )
  Model1.Conf$addMonitors( c("sigma"), print=F )
  Model1.MCMC      <- buildMCMC    ( Model1.Conf )
  Model1.MCMC.Comp <- compileNimble( Model1.MCMC )
  ntheta           <- 1e4 ; nburnin <- 1e3 ; niter <- ntheta + nburnin

  set.seed(1)
  theta <- runMCMC( Model1.MCMC.Comp, nburnin=nburnin, niter=niter, prog=F )

  alpha.NL <- theta[,"alpha"]
  beta.NL  <- theta[,"beta" ]
  sigma.NL <- theta[,"sigma"]

  par( mfrow=c(1,3), mar=c(3,3,2,0) )
  plot(alpha.NL, type="l", ylab="", main="alpha", lwd=0.1 )
  plot(beta.NL , type="l", ylab="", main="beta" , lwd=0.1 )
  plot(sigma.NL, type="l", ylab="", main="sigma", lwd=0.1 )

  par( mfrow=c(1,3), mar=c(5,5,1,1) )
  hist( alpha.NL, xlab="alpha", main="" )
  hist( beta.NL , xlab="beta" , main="" )
  hist( sigma.NL, xlab="sigma", main="" )

### PRA

  n_unc     <- 1e3
  itheta    <- sample( 1:ntheta, n_unc, replace=(ntheta<n_unc) )
  thr       <- seq( 0.1, 2.9, by=0.05 )    ;   n_thr <- length(thr)
  Ez_Hj     <- Ez_notHj  <-   Rj   <-   Vj           <- numeric(n_unc)
  Ez_H      <- Ez_notH   <-   R    <-   V    <-   pH <- numeric(n_thr)
                            s_R    <- s_V    <- s_pH <- numeric(n_thr)
  LCIR      <- UCIR      <- LCIV   <- UCIV           <- numeric(n_thr)
  LCIz_notH <- UCIz_notH <- LCIz_H <- UCIz_H         <- numeric(n_thr)
  alpha <- theta[,1] ; beta <- theta[,2] ; sigma <- theta[,3]
  for(i in 1:n_thr) {
    i_H   <- which( x < thr[i] ) ; n_H     <- length(i_H)
    pH[i] <- n_H / n             ; s_pH[i] <- sqrt( pH[i]*(1-pH[i]) / n )
    for(j in 1:n_unc) {
      zj       <- alpha[itheta[j]] + beta[itheta[j]] * exp(-x) +
                  rnorm( n, 0, sigma[itheta[j]] )
      Ez_Hj[j] <- mean( zj[i_H] ) ; Ez_notHj[j] <- mean( zj[-i_H] )
      Vj[j]    <- Ez_notHj[j] - Ez_Hj[j]
      Rj[j]    <- pH[i] * Vj[j] }
    R[i]         <- mean( Rj )     ; V[i]         <- mean( Vj )
    s_R[i]       <- sd  ( Rj )     ; s_V[i]       <- sd  ( Vj )
    Ez_H[i]      <- mean( Ez_Hj )  ; Ez_notH[i]   <- mean( Ez_notHj ) 
    qu           <- function( z, q=0.025 ){ quantile( z, q, na.rm=T ) }
    LCIz_notH[i] <- qu( Ez_notHj ) ; UCIz_notH[i] <- qu( Ez_notHj, 0.975 )
    LCIz_H[i]    <- qu( Ez_Hj    ) ; UCIz_H[i]    <- qu( Ez_Hj   , 0.975 )
    LCIV[i]      <- qu( Vj      )  ; UCIV[i]      <- qu( Vj      , 0.975 )
    LCIR[i]      <- qu( Rj      )  ; UCIR[i]      <- qu( Rj      , 0.975 )
  }

  par    ( mfrow=c(1,2), mar=c(5,5,1,1) )

  plot   ( thr, Ez_H, type="l",
           ylim=range(c(LCIz_notH,LCIz_H,UCIz_notH,UCIz_H)),
           xlab="Threshold", ylab="Expectation of z")
  polygon( c(thr,thr[n_thr:1]), c(LCIz_H,UCIz_H[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIz_notH,UCIz_notH[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  abline ( h=m[2], col="blue", lty=3 )
  lines  ( thr, Ez_H, col="red") ; lines( thr, Ez_notH, col="black" )
  legend ( "bottomright", legend=c("E[z|¬H]","E[z]","E[z|H]"),
           col=c("black","blue","red"), lty=c(1,3,1), cex=0.75 )

  plot   ( thr, V, col="black", type="l",
           ylim=range(c(LCIV,UCIV,LCIR,UCIR)),
           xlab="Threshold", ylab="V, R" )
  polygon( c(thr,thr[n_thr:1]), c(LCIV,UCIV[n_thr:1]),
           col=alpha("grey",0.6), border=NA )
  polygon( c(thr,thr[n_thr:1]), c(LCIR,UCIR[n_thr:1]),
           col=alpha("pink",0.6), border=NA )
  lines  ( thr, V, col="black" ) ; lines( thr, R, col="red" )
  legend ( "bottomright", legend=c("V","R"), col=c("black","red"), lty=1, cex=0.75 )
