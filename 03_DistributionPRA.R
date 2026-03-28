# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 3: Distribution-Based PRA

  library(copula)
  library(MCMCpack)
  library(mvtnorm)
  library(terra)
  library(VineCopula)

  load( "data/l_xz.L.RData" ) ; n_d <- length(l_xz.L)
  xz <- l_xz.L[[1]]           ; x   <- xz[,1]         ; z <- xz[,2]
  n  <- dim( xz )[1] 

## Basic Equations

  ### Verifying the analytical solutions against numerical integration
  ### (not shown in the book)
  n    <- 1e3           ;  set.seed(1)
  x    <- rbeta(n,4,2)  ;  z <- 3 + 5 * x^5 + rnorm(n,0,1)
  EzGx <- function(x) {3 + 5*x^5}
  plot(x,z)
  EzGxBelowThr.nu <- function(x,z,thr) {mean(z[x<thr])}
  EzGxBelowThr.an <- function(thr,dx=0.001) {
    intgrl <- 0  ;  ni     <- thr/dx
    for(i in 1:ni) {
      xi     <- dx * (i-0.5)  ;  di <- dbeta(xi,4,2) * EzGx(xi) * dx
      intgrl <- intgrl + di}
    return( intgrl / pbeta(thr,4,2) )
  }
  thr <- 0.4  ; c( EzGxBelowThr.nu(x,z,thr) , EzGxBelowThr.an(thr) )
  thr <- 0.7  ; c( EzGxBelowThr.nu(x,z,thr) , EzGxBelowThr.an(thr) )
  thr <- 1    ; c( EzGxBelowThr.nu(x,z,thr) , EzGxBelowThr.an(thr) )

## Fully Specified Gaussian $p[x,z]$

  mx     <- mz <- 0 ; Sx <- Sz <- 1 ; r <- 0.5
  m      <- c(mx,mz)
  S      <- diag( c(Sx,Sz) ) ; S[1,2] <- S[2,1] <- r * sqrt(Sx * Sz)
  px     <- function( x , m.=mx, S.=Sx  ) { dnorm  ( x , m., sqrt(S.) ) }
  pz     <- function( z , m.=mz, S.=Sz  ) { dnorm  ( z , m., sqrt(S.) ) }
  pxz    <- function( xz, m.=m , S.=S   ) { dmvnorm( xz, m.,      S.  ) }
  pa_b   <- function( a, b, ma, mb, Sa, Sb, r ) { dnorm( a,
    mean = ma + r*(b-mb)*sqrt(Sa/Sb), sd = sqrt( Sa*(1-r^2) ) ) }
  pz_x   <- function( z, x ) { pa_b( z, x, mz, mx, Sz, Sx, r ) }
  px_z   <- function( x, z ) { pa_b( x, z, mx, mz, Sx, Sz, r ) }

  n1      <- 41
  x.seq   <- seq( mx-2, mx+2, length.out=n1 )
  z.seq   <- seq( mz-2, mz+2, length.out=n1 )
  xz.grid <- as.matrix( expand.grid(x.seq,z.seq) )
  p.grid  <- pxz( xz.grid )
  r.xzp   <- rast( cbind( xz.grid, p.grid ), type="xyz" )
  
  par( mfrow=c(2,2), mar=c(4,2,2,1) )
  plot( x.seq, px(x.seq), type="l", main="p[x]", xlab="x", ylab="" )
  plot( z.seq, pz(z.seq), type="l", main="p[z]", xlab="z", ylab="" )
  plot( r.xzp,
        main="p[x,z]", breaks=seq(0,max(p.grid),length.out=n1),
        col=rev(gray.colors(n1,start=0,end=1)),
        axes=T, box=T, legend=F, xlab="x", ylab="z" )
  plot  ( z.seq, pz_x(z.seq,mx-1), type="l",
        main=paste0("p[z|x]"), xlab="z", ylab="" )
  points( z.seq, pz_x(z.seq,mx  ), type="l", col="red" )
  points( z.seq, pz_x(z.seq,mx+1), type="l", col="blue" )
  legend( "bottom", legend=c("x=-1","x=0","x=1"),
          col=c("black","red","blue"), lty=1, cex=0.75 )

### Hazard Probability

  thr.seq <- seq( mx-2, mx+2, length.out=101 )
  par( mfrow=c(1,2), mar=c(4,2,2,1) )
  plot( thr.seq, pnorm(thr.seq,0,1), type="l",
        main="p[H] = p[x<thr]", xlab="thr", ylab="")

  pz_xbelow <- function(thr,z, mz.=mz, mx.=mx, Sz.=Sz, Sx.=Sx, r.=r) {
    mx_z  <- mx. + r. * (z-mz.) * sqrt(Sx./Sz.)
    Sx_z  <- Sx. * (1-r.^2)
    return( pz(z) * pnorm(thr,mx_z,sqrt(Sx_z)) /
                    pnorm(thr,mx. ,sqrt(Sx. ))   ) }
  plot  ( z.seq, pz_xbelow(mx-2,z.seq), type="l", col="black",
          xlab="z", ylab="",
          main=paste0("p[z|x < thr]") )
  points( z.seq, pz_xbelow(mx  ,z.seq), type="l", col="red"  )
  points( z.seq, pz_xbelow(mx+2,z.seq), type="l", col="blue" )
  legend( "topright", legend=c("thr=-2","thr=0","thr=2"),
          col=c("black","red","blue"), lty=1, cex=0.75 )

### Vulnerability and Risk

  Ez_Gauss <- function( m.=m, S.=S, thr.=thr ) {
   mx     <- m.[1]   ; mz <- m.[2]
   Sx     <- S.[1,1] ; Sz <- S.[2,2]   ; Sxz    <- S.[1,2]
   pthr   <- dnorm(thr., mx, sqrt(Sx)) ; Fthr   <- pnorm(thr., mx, sqrt(Sx))
   Ez_xlo <- mz - Sxz * pthr / Fthr    ; Ez_xhi <- mz + Sxz * pthr / (1-Fthr)
   result <- c( Ez_xlo, Ez_xhi ) ; names ( result ) <- c( "Ez_xlo", "Ez_xhi" )
   return( result ) }

  thr.seq       <- seq( mx-2, mx+2, length.out=101 )
  Ez_xlo.seq <- sapply( thr.seq, function(thr){Ez_Gauss(m,S,thr)["Ez_xlo"]} )
  Ez_xhi.seq <- sapply( thr.seq, function(thr){Ez_Gauss(m,S,thr)["Ez_xhi"]} )

  par( mfrow=c(1,2), mar=c(4,2,2,1) )
  plot  ( thr.seq, Ez_xhi.seq, col="black", type="l",
          xlab="thr", ylab="", ylim=range(Ez_xhi.seq,Ez_xlo.seq) )
  points( thr.seq, rep(mz,101)  , col="blue" , type="l", lty=2 )
  points( thr.seq, Ez_xlo.seq, col="red"  , type="l" )
  legend( "bottomright", legend=c("E[z|x>=thr]","E[z]","E[z|x<thr]"),
          col=c("black","blue","red"),
          lty=c(1,2,1), cex=0.75 )

  V.seq <- Ez_xhi.seq - Ez_xlo.seq
  R.seq <- Ez_xhi.seq - mz
  plot  ( thr.seq, V.seq, col="black", type="l",
          ylim=c(0,max(V.seq)), xlab="thr", ylab="" )
  points( thr.seq, R.seq, col="red"  , type="l" )
  legend( "bottomright", legend=c("V","R"), col=c("black","red"),
          lty=1, cex=0.75 )

### Simple Equations for Bivariate Gaussian PRA
  
  PRA0_Gauss <- function( m.=m, S.=S, thr.=thr ) {
    mx <- m.[1] ; sx <- sqrt(S.[1,1]) ; Sxz <- S.[1,2]
    pH <- pnorm( thr., mx, sx )
    V  <- Sxz * dnorm(thr., mx, sx) / (pH * (1-pH))
    R  <- pH * V
    return( c( pH=pH, V=V, R=R ) ) }

## Gaussian $p[x,z]$ with Unknown Hyperparameters

  xz <- l_xz.L[[1]] ; n <- dim(xz)[1] ; S <- cov(xz)

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

   PRA( x, z, thr )
   PRA_Gauss( m.=m, S.=S, n.=n, thr.=thr )

  set.seed(1)

  thr <- -1
       round( PRA( x, z, thr ), 7 )
  cat( round( PRA_Gauss( m.=m, S.=S, n.=n, thr.=thr ), 7 ) )

### Numerical Check of UQ

  PRA.tbl  <- t( sapply( 1:n_d, function(d){
    PRA_Gauss( colMeans(l_xz.L[[d]]), cov(l_xz.L[[d]]), n, -1 ) } ) )

  # "True" mean and sd of the PRA-results over the n_d samples
  m_pH <- mean(PRA.tbl[,"pH"]) ; s_pH <- sd(PRA.tbl[,"pH"]) 
  m_V  <- mean(PRA.tbl[,"V" ]) ; s_V  <- sd(PRA.tbl[,"V" ]) 
  m_R  <- mean(PRA.tbl[,"R" ]) ; s_R  <- sd(PRA.tbl[,"R" ])

  par( mfrow=c(2,3), mar=c(2,2,2,2) )

  range.pH <- range( 0, PRA.tbl [,"pH"] )
  range.V  <- range( 0, PRA.tbl [,"V" ] )
  range.R  <- range( 0, PRA.tbl [,"R" ] )
  hist( PRA.tbl[,"pH"], xlim=range.pH, main="p[H]", xlab="", ylab="" )
    abline( v=m_pH, col="red", lwd=1 )
  hist( PRA.tbl[,"V" ], xlim=range.V , main="V"   , xlab="", ylab="" )
    abline( v=m_V , col="red", lwd=1 )
  hist( PRA.tbl[,"R" ], xlim=range.R , main="R"   , xlab="", ylab="" )
    abline( v=m_R , col="red", lwd=1 )
    
  range.s_pH <- range( 0, PRA.tbl[,"s_pH"] )
  range.s_V  <- range( 0, PRA.tbl[,"s_V" ] )
  range.s_R  <- range( 0, PRA.tbl[,"s_R" ] )
  hist( PRA.tbl[,"s_pH"], xlim=range.s_pH, main="sd( p[H] )", xlab="", ylab="" )
    abline( v=s_pH, col="red", lwd=1 )
  hist( PRA.tbl[,"s_V" ], xlim=range.s_V , main="sd( V )"   , xlab="", ylab="" )
    abline( v=s_V , col="red", lwd=1 )
  hist( PRA.tbl[,"s_R" ], xlim=range.s_R , main="sd( R )"   , xlab="", ylab="" )
    abline( v=s_R , col="red", lwd=1 )

### Bayesian and Non-Bayesian Alternatives

  PRA_Gauss.s <- function( m.=m, S.=S, n.=n, thr.=thr ) {
    pH <- V <- R <- rep( NA, 1e3 )
    for(j in 1:1e3){
      m     <- rmvnorm( 1, m., S./n. )
      S     <- rWishart( 1, n.-1, S./(n.-1) )[,,1]
      PRA   <- PRA0_Gauss( m, S, thr. )
      pH[j] <- PRA["pH"] ; V[j] <- PRA["V"] ; R[j] <- PRA["R"]
    }
    return( c( pH=mean(pH),   V=mean(V),   R=mean(R),
             s_pH=sd  (pH), s_V=sd  (V), s_R=sd  (R) ) )
  }

  PRA.tbl.s  <- t( sapply( 1:n_d, function(d){
    PRA_Gauss.s ( colMeans(l_xz.L[[d]]), cov(l_xz.L[[d]]), n, -1 ) } ) )

  m_pH.s <- mean(PRA.tbl.s[,"pH"]) ; s_pH.s <- sd(PRA.tbl.s[,"pH"])
  m_V.s  <- mean(PRA.tbl.s[,"V" ]) ; s_V.s  <- sd(PRA.tbl.s[,"V" ])
  m_R.s  <- mean(PRA.tbl.s[,"R" ]) ; s_R.s  <- sd(PRA.tbl.s[,"R" ])

  par( mfrow=c(2,3), mar=c(2,2,2,2) )

  range.pH.s <- range( 0, PRA.tbl.s[,"pH"] )
  range.V.s  <- range( 0, PRA.tbl.s[,"V" ] )
  range.R.s  <- range( 0, PRA.tbl.s[,"R" ] )
  hist( PRA.tbl.s[,"pH"], xlim=range.pH.s, main="p[H]", xlab="", ylab="" )
    abline( v=m_pH.s, col="red", lwd=1 )
  hist( PRA.tbl.s[,"V" ], xlim=range.V.s , main="V"   , xlab="", ylab="" )
    abline( v=m_V.s , col="red", lwd=1 )
  hist( PRA.tbl.s[,"R" ], xlim=range.R.s , main="R"   , xlab="", ylab="" )
    abline( v=m_R.s , col="red", lwd=1 )
    
  range.s_pH.s <- range( 0, PRA.tbl.s[,"s_pH"] )
  range.s_V.s  <- range( 0, PRA.tbl.s[,"s_V" ] )
  range.s_R.s  <- range( 0, PRA.tbl.s[,"s_R" ] )
  hist( PRA.tbl.s[,"s_pH"], xlim=range.s_pH.s, main="sd( p[H] )", xlab="", ylab="" )
    abline( v=s_pH.s, col="red", lwd=1 )
  hist( PRA.tbl.s[,"s_V" ], xlim=range.s_V.s , main="sd( V )"   , xlab="", ylab="" )
    abline( v=s_V.s , col="red", lwd=1 )
  hist( PRA.tbl.s[,"s_R" ], xlim=range.s_R.s , main="sd( R )"   , xlab="", ylab="" )
    abline( v=s_R.s , col="red", lwd=1 )

## Copulas for Distribution-Based PRA with Known Marginal Distributions
### Sampling from Copulas and Carrying out PRA

  cpN    <- normalCopula( param=0.5, dim=2 )
  mvN.NN <- mvdc( cpN, margins = c("norm", "norm"),
                  paramMargins = list( list( mean=0, sd=1 ),
                                       list( mean=2, sd=1 ) ) )

  n             <- 1e3
  sample.cpN    <- rCopula( n, cpN )
           Fx.N <- sample.cpN[,1]    ; Fz.N   <- sample.cpN[,2]
  sample.mvN.NN <- rMvdc  ( n, mvN.NN )
         x.N.NN <- sample.mvN.NN[,1] ; z.N.NN <- sample.mvN.NN[,2]

  par( mfrow=c(2,2), mar=c(4,4,3,0) )
  
  plot( Fx.N  , Fz.N  , main="Sample from\nGaussian copula", xlab="Fx", ylab="Fz" )
  plot( x.N.NN, z.N.NN, main="Sample from\np[x,z]", xlab="x" , ylab="z"  )
  hist( z.N.NN, main="", xlab="z", ylab="")
  hist( x.N.NN, main="", xlab="x", ylab="")

  # Verifying that the samples from copula and from p[x,z] are consistent:
  # plot( pnorm( x), pnorm( z,2),
  #       main=paste0( "r=", signif( cor(pnorm( x),pnorm( z,2)), 2 ) ),
  #       xlab="F(x)", ylab="F(z)" )
  # plot( qnorm(Fx), qnorm(Fz,2),
  #       main=paste0( "r=", signif( cor(qnorm(Fx),qnorm(Fz,2)), 2 ) ),
  #       xlab="d(Fx)" , ylab="d(Fz)"  )

  mvN.NG        <- mvdc( cpN, margins = c("norm", "gamma"),
                         paramMargins = list( list( mean =0, sd   =1 ),
                                        list( shape=4, scale=2 ) ) )
  sample.mvN.NG <- rMvdc( n, mvN.NG )
         x.N.NG <- sample.mvN.NG[,1] ; z.N.NG <- sample.mvN.NG[,2]

  par( mfrow=c(2,2), mar=c(4,4,3,0) )
  
  plot( Fx.N  , Fz.N  , main="Sample from\nGaussian copula", xlab="Fx", ylab="Fz" )
  plot( x.N.NG, z.N.NG, main="Sample from\np[x,z]", xlab="x" , ylab="z"  )
  hist( z.N.NG, main="", xlab="z", ylab="")
  hist( x.N.NG, main="", xlab="x", ylab="")

  signif( PRA( x.N.NN, z.N.NN, thr=-1 ), 2 )
  signif( PRA( x.N.NG, z.N.NG, thr=-1 ), 2 )

  cpt <- tCopula( param=0.5, df=1 )

  sample.cpt <- rCopula( n, cpt )
        Fx.t <- sample.cpt[,1] ; Fz.t   <- sample.cpt[,2]
  mvt.NG <- mvdc( cpt, margins = c("norm", "gamma"),
                  paramMargins = list( list( mean =0, sd   =1 ),
                                       list( shape=4, scale=2 ) ) )
  sample.mvt.NG <- rMvdc( n, mvt.NG )
         x.t.NG <- sample.mvt.NG[,1] ; z.t.NG <- sample.mvt.NG[,2]

  par( mfrow=c(2,2), mar=c(4,4,3,0) )
  plot( Fx.t  , Fz.t  , main="Sample from\nt-copula", xlab="Fx", ylab="Fz" )
  plot( x.t.NG, z.t.NG, main="Sample from\np[x,z]", xlab="x" , ylab="z"  )
  hist( z.t.NG, main="", xlab="z", ylab="")
  hist( x.t.NG, main="", xlab="x", ylab="")

  signif( PRA( x.t.NG, z.t.NG, thr=-1 ), 2 )
  
### Copula Selection
  
fitC <- function(x,z){ BiCopSelect( ecdf(x)(x)*n/(n+1), ecdf(z)(z)*n/(n+1),
                                    sel="BIC" ) }
fitC( x.N.NN, z.N.NN )
fitC( x.N.NG, z.N.NG )
fitC( x.t.NG, z.t.NG )

s.cpN1 <- rCopula( n, normalCopula (param=0.1) )
s.cpN2 <- rCopula( n, normalCopula (param=0.5) )
s.cpN3 <- rCopula( n, normalCopula (param=0.9) )

s.cpt1 <- rCopula( n, tCopula      (param=0.5, df= 1 ) )
s.cpt2 <- rCopula( n, tCopula      (param=0.9, df= 1) )
s.cpt3 <- rCopula( n, tCopula      (param=0.9, df=99) )

s.cpC <- rCopula( n, claytonCopula(param=10) )
s.cpF <- rCopula( n, frankCopula  (param=10) )
s.cpG <- rCopula( n, gumbelCopula (param=10) )

par( mfrow=c(3,3), mar=c(4,4,3,0) )

plot( s.cpN1[,1], s.cpN1[,2], main="Gaussian copula\nrho = 0.1",
      xlab="Fx", ylab="Fz" )
plot( s.cpN2[,1], s.cpN2[,2], main="Gaussian copula\nrho = 0.5",
      xlab="Fx", ylab="Fz" )
plot( s.cpN3[,1], s.cpN3[,2], main="Gaussian copula\nrho = 0.9",
      xlab="Fx", ylab="Fz" )

plot( s.cpt1[,1], s.cpt1[,2], main="t-copula\nrho=0.5, df=1",
      xlab="Fx", ylab="Fz" )
plot( s.cpt2[,1], s.cpt2[,2], main="t-copula\nrho=0.9, df=1",
      xlab="Fx", ylab="Fz" )
plot( s.cpt3[,1], s.cpt3[,2], main="t-copula\nrho=0.9, df=99",
      xlab="Fx", ylab="Fz" )

plot( s.cpC[,1], s.cpC[,2], main="Clayton copula\ntheta = 10",
      xlab="Fx", ylab="Fz" )
plot( s.cpF[,1], s.cpF[,2], main="Frank copula\ntheta = 10",
      xlab="Fx", ylab="Fz" )
plot( s.cpG[,1], s.cpG[,2], main="Gumbel copula\ntheta = 10",
      xlab="Fx", ylab="Fz" ) 
