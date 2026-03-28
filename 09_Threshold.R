# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 9: On Threshold Choice

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

  PRAm <- function( x, z, thr=-1:1 ) {
    n   <- length(x) ; n_thr <- length(thr)
    H   <- vector("list",n_thr)
    n_H <- pH <- V <- R <- s_pH <- s_V <- s_R <- rep(NA,n_thr)
    H[[1]]  <- which( x < thr[1] ) ; n_H[1] <- length(H[[1]])
    for(i in 2:n_thr) { H[[i]] <- which( thr[i-1] <= x & x < thr[i])
                       n_H[i] <- length(H[[i]]) } ; n_notH <- n - sum(n_H)
    H.all   <- which( x < thr[n_thr] )
    pH      <- n_H / n             ; s_pH      <- sqrt( pH*(1-pH) / n )
    Ez_notH <- mean( z[-H.all] )   ; s_Ez_notH <- sqrt( var(z[-H.all] ) / n_notH )   
    for(i in 1:n_thr) {
      Ez_Hi <- mean( z[ H[[i]] ] ) ; s_Ez_Hi <- sqrt( var(z[ H[[i]]]) /  n_H[i] )
      V[i]  <- Ez_notH - Ez_Hi     ; s_V[i]  <- sqrt( s_Ez_notH^2 + s_Ez_Hi^2 ) }
    R       <- pH * V
    s_R     <- sqrt( s_pH^2 * s_V^2 + s_pH^2 * V^2 + pH^2 * s_V^2 )
    return( cbind( thr, pH, V, R, s_pH, s_V, s_R ) )
  }

## Conditions for V Being Constant in Single-Threshold PRA
  
  set.seed(1)

  n   <- 5e3
  x_U <- runif(n)     ; z_U <- 0.5 + x_U/2               + rnorm(n,0,0.05)
  x_E <- rexp(n)      ; z_E <- 1   - exp(-x_E)/2         + rnorm(n,0,0.05)
  x_B <- rbeta(n,5,1) ; z_B <- 0.5 + pbeta(x_B,5,1)/2    + rnorm(n,0,0.05)
  x_t <- rt(n,1,30)   ; z_t <- 0.5 + pt(x_t,1,30)/2      + rnorm(n,0,0.05)
  x_G <- rnorm(n)     ; z_G <- 0.5 + pnorm(x_G)/2        + rnorm(n,0,0.05)
  x_L <- rnorm(n)     ; z_L <- 0.5 + 0.5/(1+exp(-2*x_L)) + rnorm(n,0,0.05)
  
  thr.seq_U <- quantile( x_U, (1:19)/20 )
  PRA.seq_U <- t( sapply( thr.seq_U, function(t){PRA(x_U,z_U,t)} ) )

  thr.seq_E <- quantile( x_E, (1:19)/20 )
  thr.seq_B <- quantile( x_B, (1:19)/20 )
  thr.seq_t <- quantile( x_t, (1:19)/20 )
  thr.seq_G <- quantile( x_G, (1:19)/20 )
  thr.seq_L <- quantile( x_L, (1:19)/20 )
  
  PRA.seq_E <- t( sapply( thr.seq_E, function(t){PRA(x_E,z_E,t)} ) )
  PRA.seq_B <- t( sapply( thr.seq_B, function(t){PRA(x_B,z_B,t)} ) )
  PRA.seq_t <- t( sapply( thr.seq_t, function(t){PRA(x_t,z_t,t)} ) )
  PRA.seq_G <- t( sapply( thr.seq_G, function(t){PRA(x_G,z_G,t)} ) )
  PRA.seq_L <- t( sapply( thr.seq_L, function(t){PRA(x_L,z_L,t)} ) )

  par(mfrow=c(2,3), mar=c(4,2,2,1))

  plot( x_U, z_U, xlab="", ylim=range(z_U,PRA.seq_U),
        pch=".", main="Uniform" )
  points( thr.seq_U, PRA.seq_U[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_U, PRA.seq_U[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_U, PRA.seq_U[,"R"] , type="b", lwd=3, col="red" )

  plot( x_E, z_E, xlab="", xlim=range(thr.seq_E), ylim=range(z_E,PRA.seq_E),
        pch=".", main="Exponential" )
  points( thr.seq_E, PRA.seq_E[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_E, PRA.seq_E[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_E, PRA.seq_E[,"R"] , type="b", lwd=3, col="red" )

  plot( x_B, z_B,  xlab="", ylim=range(z_B,PRA.seq_B),
          pch=".", main="Beta" )
  points( thr.seq_B, PRA.seq_B[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_B, PRA.seq_B[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_B, PRA.seq_B[,"R"] , type="b", lwd=3, col="red" )

  plot( x_t, z_t,  xlab="x or thr", xlim=range(thr.seq_t), ylim=range(z_t,PRA.seq_t),
        pch=".", main="t" )
  points( thr.seq_t, PRA.seq_t[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_t, PRA.seq_t[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_t, PRA.seq_t[,"R"] , type="b", lwd=3, col="red" )

  plot( x_G, z_G,  xlab="x or thr", ylim=range(z_G,PRA.seq_G),
        pch=".", main="Gaussian" )
  points( thr.seq_G, PRA.seq_G[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_G, PRA.seq_G[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_G, PRA.seq_G[,"R"] , type="b", lwd=3, col="red" )

  plot( x_L, z_L,  xlab="x or thr", ylim=range(z_L,PRA.seq_L),
        pch=".", main="Gaussian-Logistic" )
  points( thr.seq_L, PRA.seq_L[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_L, PRA.seq_L[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_L, PRA.seq_L[,"R"] , type="b", lwd=3, col="red" )
  legend ( "topleft", legend=c("z","pH","V","R"),
           col=c("black","blue","green","red"),
           pch=c(20,NA,NA,NA), lty=c(NA,1,1,1), cex=0.75 )
  
## Single- vs. Multi-Threshold PRA

  thr.seqU_U  <- seq( quantile(x_U,1/10), quantile(x_U,9/10), length.out=9 )
  PRAm.seqU_U <- PRAm( x_U, z_U, thr.seqU_U )

  thr.seqU_E  <- seq( quantile(x_E,1/10), quantile(x_E,9/10), length.out=9 )
  thr.seqU_B  <- seq( quantile(x_B,1/10), quantile(x_B,9/10), length.out=9 )
  thr.seqU_t  <- seq( quantile(x_t,1/10), quantile(x_t,9/10), length.out=9 )
  thr.seqU_G  <- seq( quantile(x_G,1/10), quantile(x_G,9/10), length.out=9 )
  thr.seqU_L  <- seq( quantile(x_L,1/10), quantile(x_L,9/10), length.out=9 )
  
  PRAm.seqU_E <- PRAm( x_E, z_E, thr.seqU_E )
  PRAm.seqU_B <- PRAm( x_B, z_B, thr.seqU_B )
  PRAm.seqU_t <- PRAm( x_t, z_t, thr.seqU_t )
  PRAm.seqU_G <- PRAm( x_G, z_G, thr.seqU_G )
  PRAm.seqU_L <- PRAm( x_L, z_L, thr.seqU_L )

  par(mfrow=c(2,3), mar=c(4,2,2,1))

  plot( x_U, z_U, xlab="", ylim=c(0,1.2),
        pch=".", main="Uniform" )
  points( thr.seqU_U, 2*PRAm.seqU_U[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seqU_U, 2*PRAm.seqU_U[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seqU_U, 2*PRAm.seqU_U[,"R"] , type="b", lwd=3, col="red" )

  plot( x_E, z_E, xlab="", xlim=range(thr.seq_E), ylim=c(0,1.2),
        pch=".", main="Exponential" )
  points( thr.seqU_E, 2*PRAm.seqU_E[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seqU_E, 2*PRAm.seqU_E[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seqU_E, 2*PRAm.seqU_E[,"R"] , type="b", lwd=3, col="red" )

  plot( x_B, z_B,  xlab="", ylim=c(0,1.2), pch=".", main="Beta" )
  points( thr.seqU_B, 2*PRAm.seqU_B[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seqU_B, 2*PRAm.seqU_B[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seqU_B, 2*PRAm.seqU_B[,"R"] , type="b", lwd=3, col="red" )

  plot( x_t, z_t,  xlab="x", xlim=range(thr.seq_t), ylim=c(0,1.2),
        pch=".", main="t" )
  points( thr.seqU_t, 2*PRAm.seqU_t[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seqU_t, 2*PRAm.seqU_t[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seqU_t, 2*PRAm.seqU_t[,"R"] , type="b", lwd=3, col="red" )

  plot( x_G, z_G,  xlab="x", ylim=c(0,1.2), pch=".", main="Gaussian" )
  points( thr.seqU_G, 2*PRAm.seqU_G[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seqU_G, 2*PRAm.seqU_G[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seqU_G, 2*PRAm.seqU_G[,"R"] , type="b", lwd=3, col="red" )

  plot( x_L, z_L,  xlab="x", ylim=c(0,1.2), pch=".", main="Gaussian-Logistic" )
  points( thr.seqU_L, 2*PRAm.seqU_L[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seqU_L, 2*PRAm.seqU_L[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seqU_L, 2*PRAm.seqU_L[,"R"] , type="b", lwd=3, col="red" )
  legend ( "topleft", legend=c("z","pH*10","V","R*10"),
           col=c("black","blue","green","red"),
           pch=c(20,NA,NA,NA), lty=c(NA,1,1,1), cex=0.75 )
