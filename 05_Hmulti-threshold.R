# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 5: Refining the Hazard I: Multi-Threshold PRA

  library(MCMCpack)
  library(mvtnorm)

  load( "data/xz_4.RData" )
  load( "data/l_xz.L.RData" )
  load( "data/l_xz.NL.RData" )

  PRA <- function( x, z, thr=0 ) {
    n       <- length(x)     ; H         <- which(x < thr) ; n_H <- length(H)
    Ez_H    <- mean( z[ H] ) ; s_Ez_H    <- sqrt( var(z[ H]) /    n_H  )
    Ez_notH <- mean( z[-H] ) ; s_Ez_notH <- sqrt( var(z[-H]) / (n-n_H) )
    pH      <- n_H / n       ; V         <- Ez_notH - Ez_H ; R  <- pH * V
    s_pH    <- sqrt( pH*(1-pH) / n )
    s_V     <- sqrt( s_Ez_H^2 + s_Ez_notH^2 )
    s_R     <- sqrt( s_pH^2*s_V^2 + s_pH^2*V^2 + pH^2*s_V^2 )
    return( c(pH=pH,V=V,R=R,s_pH=s_pH,s_V=s_V,s_R=s_R) ) }

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
             s_pH=sd  (pH), s_V=sd  (V), s_R=sd  (R) ) ) }

## Sampling-Based Multi-Threshold PRA

  PRAm <- function( x, z, thr=-1:1 ) {
    n   <- length(x) ; n_thr <- length(thr)
    H   <- vector("list",n_thr)
    n_H <- pH <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_thr)
    H[[1]]  <- which( x < thr[1] ) ; n_H[1] <- length(H[[1]])
    for(i in 2:n_thr) { H[[i]] <- which( thr[i-1] <= x & x < thr[i])
                       n_H[i] <- length(H[[i]]) } ; n_notH <- n - sum(n_H)
    H.all   <- which( x < thr[n_thr] )
    pH      <- n_H / n           ; s_pH      <- sqrt(pH*(1-pH) / n)
    Ez_notH <- mean( z[-H.all] ) ; s_Ez_notH <- sqrt(var(z[-H.all] ) / n_notH)   
    for(i in 1:n_thr) {
      Ez_Hi <- mean( z[ H[[i]] ] ) ; s_Ez_Hi <- sqrt(var(z[ H[[i]]]) / n_H[i])
      V[i]  <- Ez_notH - Ez_Hi     ; s_V[i]  <- sqrt(s_Ez_notH^2 + s_Ez_Hi^2) }
    R       <- pH * V
    s_R     <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2 )
    return( cbind( thr, pH, V, R, s_pH, s_V, s_R ) ) }

### Uncertainty in $p[H]$: Alternative Bayesian Approach

  pHm_Di <- function( x, thr=-1:1 ){
    n      <- length(x) ; n_thr <- length(thr)
    H      <- vector("list",n_thr)
    n_H    <- pH <- s_pH <- rep(NA,n_thr)
    H[[1]] <- which( x < thr[1] ) ; n_H[1] <- length(H[[1]])
    for(i in 2:n_thr) { H[[i]] <- which( thr[i-1] <= x & x < thr[i])
                       n_H[i] <- length(H[[i]]) }
    n_NotH <- n - sum(n_H)
    a0     <- rep(1,n_thr+1) ; a1 <- a0 + c(n_H,n_NotH) ; A1 <- sum(a1)
    pH     <- a1[1:n_thr] / A1
    s_pH   <- sqrt( (a1[1:n_thr]/A1) * (1 - a1[1:n_thr]/A1) / (A1+1) )
  return( cbind( pH, s_pH ) )
  }

## Distribution-Based Multi-Threshold PRA

  Ez_Gauss <- function( m.=m, S.=S, thr.=thr ) {
    mx     <- m.[1]   ; mz <- m.[2]
    Sx     <- S.[1,1] ; Sz <- S.[2,2]   ; Sxz    <- S.[1,2]
    pthr   <- dnorm(thr., mx, sqrt(Sx)) ; Fthr   <- pnorm(thr., mx, sqrt(Sx))
    Ez_xlo <- mz - Sxz * pthr / Fthr    ; Ez_xhi <- mz + Sxz * pthr / (1-Fthr)
    result <- c( Ez_xlo, Ez_xhi ) ; names ( result ) <- c( "Ez_xlo", "Ez_xhi" )
    return( result ) }

  PRAm_Gauss <- function( m.=m, S.=S, n.=n, thr.=-1:0 ) {
    n_thr <- length(thr.)
    pH    <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_thr)
    for(i in 1:n_thr){
      pHj <- Vj <- Rj <- rep( NA, 1e3 )
      for(j in 1:1e3){
        S       <- riwish( n.-1, S. * (n.-1) ) ; m <- rmvnorm( 1, m., S/n. )
        mx      <- m[1] ; sx <- sqrt(S[1,1]) ; Vxz <- S[1,2]
        Ez_NotH <- Ez_Gauss(m,S,thr.[n_thr])["Ez_xhi"]
        if(i==1){
          pHj[j] <- pnorm( thr.[i], mx, sx )
          Ez_Hi  <- Ez_Gauss(m,S,thr.[i])["Ez_xlo"]
        }else{
          pHj[j] <- pnorm( thr.[i], mx, sx ) - pnorm( thr.[i-1], mx, sx )
          Ez_Hi  <- (
            pnorm(thr.[i  ],mx,sx) * Ez_Gauss(m,S,thr.[i  ])["Ez_xlo"] -
            pnorm(thr.[i-1],mx,sx) * Ez_Gauss(m,S,thr.[i-1])["Ez_xlo"] ) / pHj[j]
        }
        Vj[j] <- Ez_NotH - Ez_Hi
        Rj[j] <- pHj[j] * Vj[j]
      }
      pH[i]   <- mean(pHj) ;   V[i] <- mean(Vj) ;   R[i] <- mean(Rj)
      s_pH[i] <- sd  (pHj) ; s_V[i] <- sd  (Vj) ; s_R[i] <- sd  (Rj)
    }
    return( cbind( thr., pH, V, R, s_pH, s_V, s_R ) )
  }

## Two Examples

  xz      <- l_xz.L[[1]] ; x <- xz[,1]       ; z <- xz[,2]
  thr_L2  <- -1:-0
  thr_L13 <- seq(-2,1,length=13)

  x_NL     <- l_xz.NL[[1]][,1]    ; z_NL <- l_xz.NL[[1]][,2]
  thr_NL14 <- seq( 0.2, 2.8, length=14 )
  thr_NL2  <- 1:2
  
  par( mfrow=c(1,2), mar=c(5,4,2,1) )
  plot(x   ,z   , pch=".", xlab="x", ylab="z", main="Linear"    )
    abline(v=thr_L13  ) ; abline(v=thr_L2 , col="red", lwd=2 )
  plot(x_NL,z_NL, pch=".", xlab="x", ylab="z", main="Nonlinear" )
    abline(v=thr_NL14 ) ; abline(v=thr_NL2, col="red", lwd=2 )

### Example I: Linear Relationship

  xz     <- l_xz.L[[1]] ; x <- xz[,1]       ; z <- xz[,2]
  n      <- length(x)   ; m <- colMeans(xz) ; S <- cov(xz)
  thr_L2 <- -1:-0
  PRAm      ( x, z,    thr_L2 )
  PRAm_Gauss( m, S, n, thr_L2 )

  xz  <- l_xz.L[[1]] ; x <- xz[,1]       ; z <- xz[,2]
  n   <- length(x)   ; m <- colMeans(xz) ; S <- cov(xz)
  thr_L2 <- -1:-0
  signif( PRAm      ( x, z, thr_L2    ), 2 )
  signif( PRAm_Gauss( m, S, n, thr_L2 ), 2 )

  thr_L1 <- 0
  PRA( x, z, thr_L1 )

  thr_L13 <- seq(-2,1,length=13) ; PRAm_L <- PRAm(x, z, thr_L13)

  pH <- PRAm_L[,"pH"] ; s_pH <- PRAm_L[,"s_pH"]
  V  <- PRAm_L[,"V" ] ; s_V  <- PRAm_L[,"s_V" ]
  R  <- PRAm_L[,"R" ] ; s_R  <- PRAm_L[,"s_R" ]
  
  par( mfrow=c(1,3), mar=c(5,3,3,0) )
  
  barpH <- barplot( PRAm_L[,"pH"], main="p[H]", names.arg=thr_L13,
                    ylab="",  ylim=c(0,max(pH+s_pH)))
    ew <- (barpH[2,1]-barpH[1,1]) / 4
    segments( barpH   , pH     , barpH   , pH+s_pH )
    segments( barpH-ew, pH+s_pH, barpH+ew, pH+s_pH )
  barV <- barplot( PRAm_L[,"V"] , main="V", names.arg=thr_L13,
                   xlab="Upper bound of interval",
                   ylab="", ylim=c(0,max(V+s_V)))
    ew <- (barV[2,1]-barV[1,1]) / 4
    segments( barV   , V    , barV   , V+s_V )
    segments( barV-ew, V+s_V, barV+ew, V+s_V )
  barR <- barplot( PRAm_L[,"R"] , main="R", names.arg=thr_L13,
                   ylab="", ylim=c(0,max(R+s_R)))
    ew <- (barR[2,1]-barR[1,1]) / 4
    segments( barR   , R    , barR   , R+s_R )
    segments( barR-ew, R+s_R, barR+ew, R+s_R )

### Example II: Nonlinear Relationship

  x_NL     <- l_xz.NL[[1]][,1] ; z_NL <- l_xz.NL[[1]][,2]
  thr_NL14 <- seq( 0.2, 2.8, length=14 )
  PRAm_NL  <- PRAm( x_NL, z_NL, thr_NL14 )

  pH <- PRAm_NL[,"pH"] ; s_pH <- PRAm_NL[,"s_pH"]
  V  <- PRAm_NL[,"V" ] ; s_V  <- PRAm_NL[,"s_V" ]
  R  <- PRAm_NL[,"R" ] ; s_R  <- PRAm_NL[,"s_R" ]
  
  par( mfrow=c(1,3), mar=c(5,3,3,0) )
  
  barpH <- barplot( pH, main="p[H]", names.arg=thr_NL14,
                    ylab="",  ylim=c(0,max(pH+s_pH)))
    ew <- (barpH[2,1]-barpH[1,1]) / 4
    segments( barpH   , pH     , barpH   , pH+s_pH )
    segments( barpH-ew, pH+s_pH, barpH+ew, pH+s_pH )
  barV <- barplot( V, main="V", names.arg=thr_NL14,
                   xlab="Upper bound of interval",
                   ylab="", ylim=c(0,max(V+s_V)))
    ew <- (barV[2,1]-barV[1,1]) / 4
    segments( barV   , V    , barV   , V+s_V )
    segments( barV-ew, V+s_V, barV+ew, V+s_V )
  barR <- barplot( R, main="R", names.arg=thr_NL14,
                   ylab="", ylim=c(0,max(R+s_R)))
    ew <- (barR[2,1]-barR[1,1]) / 4
    segments( barR   , R    , barR   , R+s_R )
    segments( barR-ew, R+s_R, barR+ew, R+s_R )

  thr_NL2  <- 1:2
  PRAm.tbl <- sapply( 1:length(l_xz.NL), function(d){ PRAm(
    l_xz.NL[[d]][,"x"], l_xz.NL[[d]][,"z"], thr=thr_NL2) }, simplify="array" )

  PRAm.tbl1 <- t(PRAm.tbl[1,,]) ; PRAm.tbl2 <- t(PRAm.tbl[2,,])
  s_pH1 <- sd(PRAm.tbl1[,"pH"]) ; s_pH2 <- sd(PRAm.tbl2[,"pH"])
  s_V1  <- sd(PRAm.tbl1[,"V" ]) ; s_V2  <- sd(PRAm.tbl2[,"V" ])
  s_R1  <- sd(PRAm.tbl1[,"R" ]) ; s_R2  <- sd(PRAm.tbl2[,"R" ])

  par( mfrow=c(2,3), mar=c(2,2,2,2) )

  range.s_pH <- range( 0, s_pH1, s_pH2, PRAm.tbl1[,"s_pH"], PRAm.tbl2[,"s_pH"] )
  range.s_V  <- range( 0, s_V1 , s_V2 , PRAm.tbl1[,"s_V" ], PRAm.tbl2[,"s_V" ] )
  range.s_R  <- range( 0, s_R1 , s_R2 , PRAm.tbl1[,"s_R" ], PRAm.tbl2[,"s_R" ] )

  hist( PRAm.tbl1[,"s_pH"], xlim=range.s_pH, main="sd( p[H1] )", xlab="", ylab="" )
    abline( v=s_pH1, col="red", lwd=1 )
  hist( PRAm.tbl1[,"s_V" ], xlim=range.s_V , main="sd( V1 )"   , xlab="", ylab="" )
    abline( v=s_V1 , col="red", lwd=1 )
  hist( PRAm.tbl1[,"s_R"] , xlim=range.s_R , main="sd( R1 )"   , xlab="", ylab="" )
    abline( v=s_R1 , col="red", lwd=1 )

  hist( PRAm.tbl2[,"s_pH"], xlim=range.s_pH, main="sd( p[H2] )", xlab="", ylab="" )
    abline( v=s_pH2, col="red", lwd=1 )
  hist( PRAm.tbl2[,"s_V" ], xlim=range.s_V , main="sd( V2 )"   , xlab="", ylab="" )
    abline( v=s_V2 , col="red", lwd=1 )
  hist( PRAm.tbl2[,"s_R" ], xlim=range.s_R , main="sd( R2 )"   , xlab="", ylab="" )
    abline( v=s_R2 , col="red", lwd=1 )

## Exercises

  PRAm  ( x_4, z_4, c(-2,0) )[,c("pH","s_pH")]
  pHm_Di( x_4     , c(-2,0) )

  PRAm  ( l_xz.L[[1]][,1], l_xz.L[[1]][,2], thr=-1:1)[,c("pH","s_pH")]
  pHm_Di( l_xz.L[[1]][,1], thr=-1:1 )
