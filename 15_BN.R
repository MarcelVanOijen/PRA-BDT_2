# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 15: Bayesian Networks for BDT

  library(mvtnorm)
  library(nimble)

  set.seed(1)

## Gaussian Bayesian Networks (GBN)
### Switching Between the Three Different Specifications

  ScondCregr <- function( W ) {
    n  <- dim(W)[1] ; S <- rep( NA, n )
    Wk <- W         ; C <- matrix( 0, nrow=n, ncol=n )
    for(k in n:2){
      ik      <- 1 : (k-1)
      S[k]    <- 1 / Wk[ k, k]
      C[k,ik] <-    -Wk[ k,ik] * S[k]
      Wk      <-     Wk[ik,ik] - as.matrix(C[k,ik]) %*% C[k,ik] / S[k] }
    S[1] <- 1 / Wk[1,1]
    return( list( Scond=S, Cregr=zapsmall(C) ) ) }

  precMatrix <- function( nodes, Scond, Cregr ){
    n <- length(nodes) ; W <- 1 / Scond[1]
    for(k in 2:n){
      ck       <- Cregr[ k, 1:(k-1) ]
      W_top    <- cbind( W * Scond[k] + ck %*% t(ck), -ck ) / Scond[k] 
      W_bottom <- cbind(                      -t(ck),   1 ) / Scond[k]
      W        <- rbind( W_top, W_bottom ) }
    rownames(W) <- colnames(W) <- nodes
    return(W) }

### Sampling from a GBN

  GaussCond <- function( mz, Sz, y ) {
    i <- 1 : ( length(mz) - length(y) )
    m  <- mz[i]   + Sz[i,-i] %*% solve(Sz[-i,-i]) %*% (y-mz[-i])
    S  <- Sz[i,i] - Sz[i,-i] %*% solve(Sz[-i,-i]) %*% Sz[-i,i]
    return( list( m=m, S=S ) ) }

## A Linear BDT Example Implemented as a GBN

  nodes <-         c( "RAIN", "IRRIG", "YC", "BENEFIT", "COST", "UTILITY")
  m     <-         c( 500   ,   0    ,  5  , 150      ,   0   , 150      )
  Scond <-         c( 1e4   ,   0.01 ,  1  , 100      , 100   , 100      )
  Cregr <- matrix( c(   0   ,   0    ,  0  ,   0      ,   0   ,   0 ,
                        0   ,   0    ,  0  ,   0      ,   0   ,   0 ,
                        0.01,   0.01 ,  0  ,   0      ,   0   ,   0 ,
                        0   ,   0    , 30  ,   0      ,   0   ,   0 ,
                        0   ,   0.5  ,  0  ,   0      ,   0   ,   0 ,
                        0   ,   0    ,  0  ,   1      ,  -1   ,   0      ),
                   byrow=T, nrow=6 )

  W <- precMatrix( nodes, Scond, Cregr )
  S <- solve(W)
  ### We can retrieve the {Scond,Cregr} GBN-specification by:
  ###   ScondCregr(W) ; ScondCregr(solve(S))
  
  ### We visualise the prior uncertainty for utility by
  ### sampling from the joint Gaussian distribution
  n       <- 1e4
  sample  <- rmvnorm( n, m, S )
  colnames( sample ) <- nodes
  sampleU <- sample[,"UTILITY"]
  sampleB <- sample[,"BENEFIT"]
  sampleC <- sample[,"COST"]
  ### Using R's 'lm'-function to confirm that our GBN treats utility
  ### as the difference between benefits and costs:
  ###   lm( sampleU ~ sampleB + sampleC )
  mU  <- round( mean(sampleU), 1 )
  sdU <- round( sd  (sampleU), 1 )
  SU  <- round( var (sampleU), 1 )

  ### New ordering of nodes to have the observed node (RAIN) at the end:
  i.alt   <- c(2:6,1)  ; nodes.alt <- nodes[i.alt]
  m.alt   <- m[i.alt]  ; S.alt     <- S    [i.alt,i.alt]
  ### Adding information about the last node (RAIN):
  GBN1    <- GaussCond( m.alt, S.alt, 600 )
  m1      <- GBN1$m ; S1 <- GBN1$S ; nodes1 <- colnames(S1)

  ### New ordering of nodes to have the RAIN and IRRIG at the end:
  i.alt2  <- c(3:6,1:2) ; nodes.alt2 <- nodes[i.alt2]
  m.alt2  <- m[i.alt2]  ; S.alt2     <- S    [i.alt2,i.alt2]

  GBN2 <- GaussCond( m.alt2, S.alt2, c(600,400) )
  m2   <- GBN2$m ; S2 <- GBN2$S ; nodes2 <- colnames(S2)

## A Linear BDT Example Implemented Using $\texttt{Nimble}$

  BDT.Code <- nimbleCode({
    RAIN     ~ dnorm( 500, sd= 100   )
    IRRIG    ~ dnorm(   0, sd=   0.1 )
    YC      <- (RAIN + IRRIG) * 0.01 + eps.YC
    BENEFIT <-  YC            * 30   + eps.BE
    COST    <-  IRRIG         * 0.5  + eps.CO
    UTILITY <-  BENEFIT - COST       + eps.UT
    eps.YC   ~ dnorm( 0, sd= 1 )
    eps.BE   ~ dnorm( 0, sd=10 )
    eps.CO   ~ dnorm( 0, sd=10 )
    eps.UT   ~ dnorm( 0, sd=10 ) } )

  BDT0.Nimble    <- nimbleModel  ( BDT.Code )
  BDT0.Comp      <- compileNimble( BDT0.Nimble )
  BDT0.Conf      <- configureMCMC( BDT0.Nimble, print=F )
  BDT0.Conf$addMonitors( c("UTILITY"), print=F )
  BDT0.MCMC      <- buildMCMC    ( BDT0.Conf )
  BDT0.MCMC.Comp <- compileNimble( BDT0.MCMC )
  
  BDT1.Data      <- list( RAIN=600 )
  BDT1.Nimble    <- nimbleModel  ( BDT.Code, data=BDT1.Data )
  BDT1.Comp      <- compileNimble( BDT1.Nimble )
  BDT1.Conf      <- configureMCMC( BDT1.Nimble, print=F )
  BDT1.Conf$addMonitors( c("UTILITY"), print=F )
  BDT1.MCMC      <- buildMCMC    ( BDT1.Conf )
  BDT1.MCMC.Comp <- compileNimble( BDT1.MCMC )
  
  BDT2.Data      <- list( RAIN=600, IRRIG=400 )
  BDT2.Nimble    <- nimbleModel  ( BDT.Code, data=BDT2.Data )
  BDT2.Comp      <- compileNimble( BDT2.Nimble )
  BDT2.Conf      <- configureMCMC( BDT2.Nimble, print=F )
  BDT2.Conf$addMonitors( c("UTILITY"), print=F )
  BDT2.MCMC      <- buildMCMC    ( BDT2.Conf )
  BDT2.MCMC.Comp <- compileNimble( BDT2.MCMC )

  nsims <- 1e4 ; nburnin <- 1e3 ; niter <- nsims+nburnin
  set.seed(1)
  
  samples0  <- runMCMC( BDT0.MCMC.Comp,
                        nburnin=nburnin, niter=niter, prog=F )
  samples0U <- samples0[ , "UTILITY" ]
  m0        <- round( mean(samples0U  ), 1 )
  S0        <- round( var (samples0U  ), 0 )

  samples1  <- runMCMC( BDT1.MCMC.Comp,
                        nburnin=nburnin, niter=niter, prog=F )
  samples1U <- samples1[ , "UTILITY" ]
  m1        <- round( mean(samples1U  ), 1 )
  S1        <- round( var (samples1U  ), 0 )

  samples2  <- runMCMC( BDT2.MCMC.Comp,
                        nburnin=nburnin, niter=niter, prog=F )
  samples2U <- samples2[ , "UTILITY"]
  m2        <- round( mean(samples2U), 1 )
  S2        <- round( var (samples2U), 0 )

  par( mfrow=c(1,3), mar=c(2,5,3,5) )

  hist( samples0[,"UTILITY"],
        main=paste0( "Prior utility 1\nm=", m0, ", S=", S0 ),
        xlab="", ylab="" )
  abline( v=m0, col="red", lwd=1 )
  
  hist( samples1[,"UTILITY"],
        main=paste0( "Posterior utility 1\nm=", m1, ", S=", S1 ),
        xlab="", ylab="" )
  abline( v=m1, col="red", lwd=1 )
  
  hist( samples2[,"UTILITY"],
        main=paste0( "Posterior utility 2\nm=", m2, ", S=", S2 ),
        xlab="", ylab="" )
  abline( v=m2, col="red", lwd=1 )

### Varying IRRIG to Identify the Value for which E[u] is Maximized

  IRRIG.seq <- seq( 0, 500, by=100 ) ; ni <- length( IRRIG.seq )
  Eu.seq    <- rep (NA, ni )
  for(i in 1:ni) {
    BDTi.Data      <- list( RAIN=600, IRRIG=IRRIG.seq[i] )
    BDTi.Nimble    <- nimbleModel  ( BDT.Code, data=BDTi.Data )
    BDTi.Comp      <- compileNimble( BDTi.Nimble )
    BDTi.Conf      <- configureMCMC( BDTi.Nimble, print=F )
    BDTi.Conf$addMonitors( c("UTILITY"), print=F )
    BDTi.MCMC      <- buildMCMC    ( BDTi.Conf )
    BDTi.MCMC.Comp <- compileNimble( BDTi.MCMC )
    samplesi       <- runMCMC( BDTi.MCMC.Comp,
                               nburnin=nburnin, niter=niter, prog=F )
    Eu.seq[i]      <- mean( samplesi[,"UTILITY"] )
  }

  par( mfrow=c(1,1), mar=c(5,5,3,3) )

  IRRIGopt <- IRRIG.seq[ which.max(Eu.seq) ]
  plot( IRRIG.seq, Eu.seq, type="b", xlab="IRRIG (mm y-1)", ylab="E [ utility ]" )
  abline( v=IRRIGopt, col="red", lty=2 )
  text( IRRIGopt, mean(Eu.seq), labels=paste("a_opt =",IRRIGopt ),
        pos=4, col="red" )

## A Nonlinear BDT Example Implemented Using $\texttt{Nimble}$

  BDT.nl.Code <- nimbleCode({
    RAIN     ~  dnorm( 500, sd=100   )
    IRRIG    ~  dnorm(   0, sd=  0.1 )
    YC       <- 50 * (1-exp(-0.001*(RAIN + IRRIG))) + eps.YC
    BENEFIT  <- YC    * 30     + eps.BE
    COST     <- IRRIG * 0.5    + eps.CO
    UTILITY  <- BENEFIT - COST + eps.UT
    eps.YC   ~  dnorm( 0, sd= 1 )
    eps.BE   ~  dnorm( 0, sd=10 )
    eps.CO   ~  dnorm( 0, sd=10 )
    eps.UT   ~  dnorm( 0, sd=10 )
  } )

  BDT.nl0.Nimble    <- nimbleModel  ( BDT.nl.Code )
  BDT.nl0.Comp      <- compileNimble( BDT.nl0.Nimble )
  BDT.nl0.Conf      <- configureMCMC( BDT.nl0.Nimble, print=F )
  BDT.nl0.Conf$addMonitors( c("UTILITY"), print=F )
  BDT.nl0.MCMC      <- buildMCMC    ( BDT.nl0.Conf )
  BDT.nl0.MCMC.Comp <- compileNimble( BDT.nl0.MCMC )
  
  BDT.nl1.Data      <- list( RAIN=600 )
  BDT.nl1.Nimble    <- nimbleModel  ( BDT.nl.Code, data=BDT.nl1.Data )
  BDT.nl1.Comp      <- compileNimble( BDT.nl1.Nimble )
  BDT.nl1.Conf      <- configureMCMC( BDT.nl1.Nimble, print=F )
  BDT.nl1.Conf$addMonitors( c("UTILITY"), print=F )
  BDT.nl1.MCMC      <- buildMCMC    ( BDT.nl1.Conf )
  BDT.nl1.MCMC.Comp <- compileNimble( BDT.nl1.MCMC )
  
  BDT.nl2.Data      <- list( RAIN=600, IRRIG=400 )
  BDT.nl2.Nimble    <- nimbleModel  ( BDT.nl.Code, data=BDT.nl2.Data )
  BDT.nl2.Comp      <- compileNimble( BDT.nl2.Nimble )
  BDT.nl2.Conf      <- configureMCMC( BDT.nl2.Nimble, print=F )
  BDT.nl2.Conf$addMonitors( c("UTILITY"), print=F )
  BDT.nl2.MCMC      <- buildMCMC    ( BDT.nl2.Conf )
  BDT.nl2.MCMC.Comp <- compileNimble( BDT.nl2.MCMC )

  nsims <- 1e4 ; nburnin <- 1e3 ; niter <- nsims+nburnin
  set.seed(1)
  
  samples.nl0  <- runMCMC( BDT.nl0.MCMC.Comp,
                        nburnin=nburnin, niter=niter, prog=F )
  samples.nl0U <- samples.nl0[ , "UTILITY" ]
  m.nl0        <- round( mean(samples.nl0U  ), 1 )
  S.nl0        <- round( var (samples.nl0U  ), 0 )

  samples.nl1  <- runMCMC( BDT.nl1.MCMC.Comp,
                        nburnin=nburnin, niter=niter, prog=F )
  samples.nl1U <- samples.nl1[ , "UTILITY" ]
  m.nl1        <- round( mean(samples.nl1U  ), 1 )
  S.nl1        <- round( var (samples.nl1U  ), 0 )

  samples.nl2  <- runMCMC( BDT.nl2.MCMC.Comp,
                          nburnin=nburnin, niter=niter, prog=F )
  samples.nl2U <- samples.nl2[ , "UTILITY"]
  m.nl2        <- round( mean(samples.nl2U), 1 )
  S.nl2        <- round( var (samples.nl2U), 0 )
  
  par( mfrow=c(1,3), mar=c(2,5,3,5) )

  hist( samples.nl0U,
        main=paste0( "Prior utility 1\nm=", m.nl0, ", S=", S.nl0 ),
        xlab="", ylab="" )
  abline( v=m.nl0, col="red", lwd=1 )
  
  hist( samples.nl1U,
        main=paste0( "Posterior utility 1\nm=", m.nl1, ", S=", S.nl1 ),
        xlab="", ylab="" )
  abline( v=m.nl1, col="red", lwd=1 )
  
  hist( samples.nl2U,
        main=paste0( "Posterior utility 2\nm=", m.nl2, ", S=", S.nl2 ),
        xlab="", ylab="" )
  abline( v=m.nl2, col="red", lwd=1 )

  IRRIG.seq <- seq( 0, 1000, by=100 ) ; ni <- length( IRRIG.seq )
  Eu.nl.seq <- rep (NA, ni )
  for(i in 1:ni) {
    BDT.nli.Data      <- list( RAIN=600, IRRIG=IRRIG.seq[i] )
    BDT.nli.Nimble    <- nimbleModel  ( BDT.nl.Code, data=BDT.nli.Data )
    BDT.nli.Comp      <- compileNimble( BDT.nli.Nimble )
    BDT.nli.Conf      <- configureMCMC( BDT.nli.Nimble, print=F )
    BDT.nli.Conf$addMonitors( c("UTILITY"), print=F )
    BDT.nli.MCMC      <- buildMCMC    ( BDT.nli.Conf )
    BDT.nli.MCMC.Comp <- compileNimble( BDT.nli.MCMC )
    samples.nli       <- runMCMC( BDT.nli.MCMC.Comp,
                               nburnin=nburnin, niter=niter, prog=F )
    Eu.nl.seq[i]      <- mean( samples.nli[,"UTILITY"] )
  }
  
  par( mfrow=c(1,1), mar=c(5,5,3,3) )

  IRRIG.nl.opt <- IRRIG.seq[ which.max(Eu.nl.seq) ]
  plot( IRRIG.seq, Eu.nl.seq, type="b",
        xlab="IRRIG (mm y-1)", ylab="E [ utility ]" )
  abline( v=IRRIG.nl.opt, col="red", lty=2 )
  text( IRRIG.nl.opt, mean(Eu.nl.seq), labels=paste("a_opt =",IRRIG.nl.opt ),
        pos=4, col="red" )

## Exercise

  m_G3 <- rep(0,3) ; S_G3 <- matrix(0.5,nrow=3,ncol=3) ; diag(S_G3) <- 1
