# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Appendix C: Solutions to Exercises

  library(terra)

  load( "data/l_xz.L.RData" )

  pH_Be <- function( x, thr=0 ) {
    n  <- length(x) ; H    <- which(x < thr) ; n_H <- length(H)
    a  <- 1 + n_H   ; b    <- 1 + n - n_H
    pH <- a / (a+b) ; s_pH <- sqrt( a*b/(a+b+1) ) / (a+b)
    return( c( pH=pH, s_pH=s_pH ) )
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
  
  set.seed(1)

  n   <- 5e3

  x_U <- runif(n)     ; z_U <- 0.5 + x_U/2               + rnorm(n,0,0.05)
  x_E <- rexp(n)      ; z_E <- 1   - exp(-x_E)/2         + rnorm(n,0,0.05)
  x_B <- rbeta(n,5,1) ; z_B <- 0.5 + pbeta(x_B,5,1)/2    + rnorm(n,0,0.05)
  x_t <- rt(n,1,30)   ; z_t <- 0.5 + pt(x_t,1,30)/2      + rnorm(n,0,0.05)
  x_G <- rnorm(n)     ; z_G <- 0.5 + pnorm(x_G)/2        + rnorm(n,0,0.05)
  x_L <- rnorm(n)     ; z_L <- 0.5 + 0.5/(1+exp(-2*x_L)) + rnorm(n,0,0.05)
  
  thr.seq_U <- quantile( x_U, (1:19)/20 )
  thr.seq_E <- quantile( x_E, (1:19)/20 )
  thr.seq_B <- quantile( x_B, (1:19)/20 )
  thr.seq_t <- quantile( x_t, (1:19)/20 )
  thr.seq_G <- quantile( x_G, (1:19)/20 )
  thr.seq_L <- quantile( x_L, (1:19)/20 )

### Chapter 2

  pH_Be( l_xz.L[[1]][,1], -1 )

### Chapter 9
  
  PRAm.seq_U <- PRAm( x_U, z_U, thr.seq_U )
  PRAm.seq_E <- PRAm( x_E, z_E, thr.seq_E )
  PRAm.seq_B <- PRAm( x_B, z_B, thr.seq_B )
  PRAm.seq_t <- PRAm( x_t, z_t, thr.seq_t )
  PRAm.seq_G <- PRAm( x_G, z_G, thr.seq_G )
  PRAm.seq_L <- PRAm( x_L, z_L, thr.seq_L )

  par(mfrow=c(2,3), mar=c(4,2,2,1))

  plot( x_U, z_U, xlab="", ylim=c(0,1.2),
        pch=".", main="Uniform" )
  points( thr.seq_U, 10*PRAm.seq_U[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_U,    PRAm.seq_U[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_U, 10*PRAm.seq_U[,"R"] , type="b", lwd=3, col="red" )

  plot( x_E, z_E, xlab="", xlim=range(thr.seq_E), ylim=c(0,1.2),
        pch=".", main="Exponential" )
  points( thr.seq_E, 10*PRAm.seq_E[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_E,    PRAm.seq_E[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_E, 10*PRAm.seq_E[,"R"] , type="b", lwd=3, col="red" )

  plot( x_B, z_B,  xlab="", ylim=c(0,1.2), pch=".", main="Beta" )
  points( thr.seq_B, 10*PRAm.seq_B[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_B,    PRAm.seq_B[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_B, 10*PRAm.seq_B[,"R"] , type="b", lwd=3, col="red" )

  plot( x_t, z_t,  xlab="x", xlim=range(thr.seq_t), ylim=c(0,1.2),
        pch=".", main="t" )
  points( thr.seq_t, 10*PRAm.seq_t[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_t,    PRAm.seq_t[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_t, 10*PRAm.seq_t[,"R"] , type="b", lwd=3, col="red" )

  plot( x_G, z_G,  xlab="x", ylim=c(0,1.2), pch=".", main="Gaussian" )
  points( thr.seq_G, 10*PRAm.seq_G[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_G,    PRAm.seq_G[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_G, 10*PRAm.seq_G[,"R"] , type="b", lwd=3, col="red" )

  plot( x_L, z_L,  xlab="x", ylim=c(0,1.2), pch=".", main="Gaussian-Logistic" )
  points( thr.seq_L, 10*PRAm.seq_L[,"pH"], type="b", lwd=3, col="blue" )
  points( thr.seq_L,    PRAm.seq_L[,"V"] , type="b", lwd=3, col="green" )
  points( thr.seq_L, 10*PRAm.seq_L[,"R"] , type="b", lwd=3, col="red" )
  legend ( "topleft", legend=c("z","p[H]*10","V","R*10"),
           col=c("black","blue","green","red"),
           pch=c(20,NA,NA,NA), lty=c(NA,1,1,1), cex=0.75 )

### Chapter 10
  
  Ri <- rnorm(1e6,0.5,0.25) * rnorm(1e6,4,2) ; c(mean(Ri), sd(Ri))
  
### Chapter 13

  set.seed(1)

  pA <- 0.9 ; pB_A <- 0.9 ; pB_a <- 0.6 ; pC_B <- 0.9 ; pC_b <- 0.2
  # Analytical solutions:
  pB   <- pA   * pB_A     + (1-pA)   * pB_a
  pC   <- pB   * pC_B     + (1-pB)   * pC_b
  pC_a <- pB_a * pC_B     + (1-pB_a) * pC_b
  pA_b <- pA   * (1-pB_A) / (1-pB)
  pc_A <- pB_A * (1-pC_B) + (1-pB_A) * (1-pC_b)
  pA_c <- pA   * pc_A     / (1-pC)
  pB_c <- pB   * (1-pC_B) / (1-pC)
  # Numerical solutions:
  n    <- 1e6
  A    <- rbinom( n, 1, pA )
  B    <- rbinom( n, 1, A*pB_A + (1-A)*pB_a )
  C    <- rbinom( n, 1, B*pC_B + (1-B)*pC_b )
  pB   <- sum(B==1) / n                ; pC   <- sum(C==1) / n
  pB_a <- sum(A==0 & B==1) / sum(A==0) ; pC_a <- sum(A==0 & C==1) / sum(A==0)
  pA_b <- sum(A==1 & B==0) / sum(B==0) ; pC_b <- sum(B==0 & C==1) / sum(B==0)
  pA_c <- sum(A==1 & C==0) / sum(C==0) ; pB_c <- sum(B==1 & C==0) / sum(C==0)
