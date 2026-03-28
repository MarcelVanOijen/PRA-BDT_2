# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Appendix B: Datasets

  library(mvtnorm)
  
### Sparse Dataset: {x_4,z_4}

  set.seed(1)
  n_4  <- 4 ; Sigma_4 <- matrix( c(1,0.5,0.5,1), nrow=2 )
  xz_4 <- rmvnorm( n_4, c(0,0), Sigma_4 ) %*% chol( Sigma_4 )
  xz_4 <- sweep( xz_4, 2, colMeans(xz_4) )
  xz_4 <- xz_4 %*% solve( chol(cov(xz_4)) ) %*% chol(Sigma_4)
  x_4  <- xz_4[,1]       ; z_4 <- xz_4[,2]
  m_4  <- colMeans(xz_4) ; S_4 <- cov(xz_4)
  
  save( x_4, z_4, file="data/xz_4.RData" )

  par( mfrow=c(1,1) )
  plot( x_4, z_4, xlab="x", ylab="z", bty="n",
        xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), pch=20 )

### Collection of Linear Datasets: l_xz.L

  set.seed(1)
  mu  <- c(0,0) ; Sigma  <- diag(1,2) ; Sigma[1,2] <- Sigma[2,1] <- 0.5
  n_d <- 1e2    ; n      <- 500       ; l_xz.L     <- vector("list",n_d)
  for(d in 1:n_d) { l_xz.L[[d]] <- rmvnorm( n, mu, Sigma ) }

  save( l_xz.L, file="data/l_xz.L.RData" )

  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  plot( l_xz.L[[1]][,1], l_xz.L[[1]][,2], main="dataset 1",
        xlab="x", ylab="z", pch="." )
  plot( l_xz.L[[2]][,1], l_xz.L[[2]][,2], main="dataset 2",
        xlab="x", ylab="z", pch="." )

### Collection of Nonlinear Datasets: l_xz.NL

  set.seed(1)
  sz  <- 0.1
  n_d <- 1e2 ; n <- 500 ; l_xz.NL <- vector("list",n_d)
  for(d in 1:n_d) {
    x <- runif( n, 0, 3 ) ; ez <- rnorm( n, 0, sz) ; z <- 1-exp(-x) + ez
    l_xz.NL[[d]] <- cbind(x,z) }

  save( l_xz.NL, file="data/l_xz.NL.RData" )

  par( mfrow=c(1,2), mar=c(4,4,1,4) )
  plot( l_xz.NL[[1]][,"x"], l_xz.NL[[1]][,"z"], main="dataset 1",
        xlab="x", ylab="z", pch="." )
  plot( l_xz.NL[[2]][,"x"], l_xz.NL[[2]][,"z"], main="dataset 2",
        xlab="x", ylab="z", pch="." )
