# Probabilistic Risk Analysis and Bayesian Decision Theory - 2nd ed.
# Marcel van Oijen & Mark Brewer (2026)

# Chapter 13: Introduction to Bayesian Decision Theory (BDT)

## Model-Based BDT
  
  u <- function( a, x, f, t, e, ka, kz ) {
    z       <- f(a,x,t) + e
    cost    <- ka * a
    benefit <- kz * z
    return( benefit - cost ) }

  fz <- function(a,x,t) { t * (1-exp(-a-x)) }
  fu <- function( a, x=1, t=1, e=0, ka=0.2, kz=1 ) {
                  u( a, x, f=fz, t, e, ka, kz) }

  par( mfrow=c(1,1), mar=c(4,2,3,0) )

  na    <- 31
  a.seq <- seq( 0, 3, length.out=na )
  u.seq <- fu(a.seq)
  imaxU <- which( u.seq == max(u.seq) ) ; amaxU <- a.seq[imaxU]
  plot( a.seq, fu(a.seq), xlab="a", ylab="", main="Utility\n(no uncertainty)",
        type="l" )
  abline( v=amaxU, col="red", lty=2 )
  text( amaxU, mean(u.seq), labels=paste("a* =",amaxU ), pos=4 )

  set.seed(1)

  np     <- 5e3
  x.smp  <- rnorm( np, 1  , 1   )
  t.smp  <- rnorm( np, 1  , 0.5 )
  e.smp  <- rnorm( np, 0  , 1   )
  ka.smp <- runif( np, 0.1, 0.3 )
  kz.smp <- runif( np, 0.5, 1.5 )
  
  Qlo    <- 0.25 ; Qhi <- 0.75

  kunc    <- c( 1/1.5, 1, 1.5 )

  x.smp1  <- 1   + (x.smp  - 1   ) * kunc[1]
  t.smp1  <- 1   + (t.smp  - 1   ) * kunc[1]
  e.smp1  <-        e.smp          * kunc[1]
  ka.smp1 <- 0.2 + (ka.smp - 0.2 ) * kunc[1]
  kz.smp1 <- 1   + (kz.smp - 1   ) * kunc[1]
  
  x.smp2  <- 1   + (x.smp  - 1   ) * kunc[2]
  t.smp2  <- 1   + (t.smp  - 1   ) * kunc[2]
  e.smp2  <-        e.smp          * kunc[2]
  ka.smp2 <- 0.2 + (ka.smp - 0.2 ) * kunc[2]
  kz.smp2 <- 1   + (kz.smp - 1   ) * kunc[2]
  
  x.smp3  <- 1   + (x.smp  - 1   ) * kunc[3]
  t.smp3  <- 1   + (t.smp  - 1   ) * kunc[3]
  e.smp3  <-        e.smp          * kunc[3]
  ka.smp3 <- 0.2 + (ka.smp - 0.2 ) * kunc[3]
  kz.smp3 <- 1   + (kz.smp - 1   ) * kunc[3]
  
  u.tbl1  <- u.tbl2 <- u.tbl3 <- NULL
  for(i in 1:np) {
    u.tbl1 <- rbind( u.tbl1,
      fu( a.seq, x.smp1[i], t.smp1[i], e.smp1[i], ka.smp1[i], kz.smp1[i] ) )
    u.tbl2 <- rbind( u.tbl2,
      fu( a.seq, x.smp2[i], t.smp2[i], e.smp2[i], ka.smp2[i], kz.smp2[i] ) )
    u.tbl3 <- rbind( u.tbl3,
      fu( a.seq, x.smp3[i], t.smp3[i], e.smp3[i], ka.smp3[i], kz.smp3[i] ) )
  }
  
  uQlo.seq1 <- sapply( 1:na, function(i){ quantile(u.tbl1[,i],probs=Qlo) } )
  uQlo.seq2 <- sapply( 1:na, function(i){ quantile(u.tbl2[,i],probs=Qlo) } )
  uQlo.seq3 <- sapply( 1:na, function(i){ quantile(u.tbl3[,i],probs=Qlo) } )
  U.seq1    <- colMeans(u.tbl1) # Naming: capital U for mean or expectation of u
  U.seq2    <- colMeans(u.tbl2)
  U.seq3    <- colMeans(u.tbl3)
  uQhi.seq1 <- sapply( 1:na, function(i){ quantile(u.tbl1[,i],probs=Qhi) } )
  uQhi.seq2 <- sapply( 1:na, function(i){ quantile(u.tbl2[,i],probs=Qhi) } )
  uQhi.seq3 <- sapply( 1:na, function(i){ quantile(u.tbl3[,i],probs=Qhi) } )

  imaxU.1 <- which( U.seq1 == max(U.seq1) ) ; amaxU.1 <- a.seq[imaxU.1]
  imaxU.2 <- which( U.seq2 == max(U.seq2) ) ; amaxU.2 <- a.seq[imaxU.2]
  imaxU.3 <- which( U.seq3 == max(U.seq3) ) ; amaxU.3 <- a.seq[imaxU.3]

  par( mfrow=c(1,3), mar=c(4,2,3,0) )

  U.range <- range( uQlo.seq1, U.seq1, uQhi.seq1,
                    uQlo.seq2, U.seq2, uQhi.seq2,
                    uQlo.seq3, U.seq3, uQhi.seq3 )
  
  plot  ( a.seq, U.seq1 , xlab="a", ylab="", type="l", 
          main="Utility\n(low uncertainty)", ylim=U.range )
  points( a.seq, uQlo.seq1, type="l", lty=2 )
  points( a.seq, uQhi.seq1, type="l", lty=2 )
  abline( v=amaxU.1, col="red", lty=2 )
  text( amaxU.1, min(U.seq1), labels=paste("a* =",amaxU.1 ), pos=4 )

  plot  ( a.seq, U.seq2 , xlab="a", ylab="", type="l", 
          main="Utility", ylim=U.range )
  points( a.seq, uQlo.seq2, type="l", lty=2 )
  points( a.seq, uQhi.seq2, type="l", lty=2 )
  abline( v=amaxU.2, col="red", lty=2 )
  text( amaxU.2, min(U.seq2), labels=paste("a* =",amaxU.2 ), pos=4 )

  plot  ( a.seq, U.seq3 , xlab="a", ylab="", type="l", 
          main="Utility\n(high uncertainty)", ylim=U.range )
  points( a.seq, uQlo.seq3, type="l", lty=2 )
  points( a.seq, uQhi.seq3, type="l", lty=2 )
  abline( v=amaxU.3, col="red", lty=2 )
  text( amaxU.3, min(U.seq3), labels=paste("a* =",amaxU.3 ), pos=2 )

################################################################################
  
  # Chapter 14: Linkages between PRA and BDT
  
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
  
  na <- length(a.seq )
  np <- length(x.smp1)
  
  thr      <- 0
  PRA.BDT1 <- sapply( 1:na, function(i){ PRA( x.smp1+a.seq[i], u.tbl1[,i], thr ) } )
  PRA.BDT2 <- sapply( 1:na, function(i){ PRA( x.smp2+a.seq[i], u.tbl2[,i], thr ) } )
  PRA.BDT3 <- sapply( 1:na, function(i){ PRA( x.smp3+a.seq[i], u.tbl3[,i], thr ) } )
  
  par( mfrow=c(1,3), mar=c(4,2,3,0) )
  pHlim <- range( c( 0, PRA.BDT1["pH",], PRA.BDT2["pH",], PRA.BDT3["pH",] ) )
  Vlim  <- range( c( 0, PRA.BDT1["V" ,], PRA.BDT2["V" ,], PRA.BDT3["V" ,] ), na.rm=T )
  
  plot  ( a.seq, PRA.BDT1["pH",], type="l", ylim=pHlim,
          main="p[H]", xlab="a", ylab="")
  points( a.seq, PRA.BDT2["pH",], type="l", lty=2 )
  points( a.seq, PRA.BDT3["pH",], type="l", lty=3 )
  legend( "topright", legend=c("High unc","", "Low unc"), lty=3:1, cex=1 )
  
  plot  ( a.seq, PRA.BDT1["V",], col="black", type="l", ylim=Vlim,
          main="V", xlab="a", ylab="" )
  points( a.seq, PRA.BDT2["V",], col="black", type="l", lty=2 )
  points( a.seq, PRA.BDT3["V",], col="black", type="l", lty=3 )
  legend( "topright", legend=c("High unc","", "Low unc"), lty=3:1, cex=1 )
  
  plot  ( a.seq, PRA.BDT1["R",], col="black", type="l", ylim=Vlim,
          main="R", xlab="a", ylab="" )
  points( a.seq, PRA.BDT2["R",], col="black", type="l", lty=2 )
  points( a.seq, PRA.BDT3["R",], col="black", type="l", lty=3 )
  legend( "topright", legend=c("High unc","", "Low unc"), lty=3:1, cex=1, col="red" )
  
  ## Utility Maximisation in BDT vs. Risk Assessment in PRA: $R_c$
  
  z.tbl1 <- z.tbl2 <- z.tbl3 <- NULL
  for(i in 1:np) {
    z.tbl1 <- rbind( z.tbl1, fz(a.seq,x.smp1[i],t.smp1[i]) + e.smp1[i] )
    z.tbl2 <- rbind( z.tbl2, fz(a.seq,x.smp2[i],t.smp2[i]) + e.smp2[i] )
    z.tbl3 <- rbind( z.tbl3, fz(a.seq,x.smp3[i],t.smp3[i]) + e.smp3[i] )
  }
  
  h.lst1 <- h.lst2 <- h.lst3 <- vector( "list", na )
  for(i in 1:na) {
    h.lst1[[i]] <- which( (a.seq[i]+x.smp1) > thr )
    h.lst2[[i]] <- which( (a.seq[i]+x.smp2) > thr )
    h.lst3[[i]] <- which( (a.seq[i]+x.smp3) > thr )
  }
  
  benefit.mnabove.seq1 <- benefit.mnabove.seq2 <- benefit.mnabove.seq3 <- rep(NA,na)
  cost.mnabove.seq1    <- cost.mnabove.seq2    <- cost.mnabove.seq3    <- rep(NA,na)
  for(i in 1:na) {
    benefit.mnabove.seq1[i] <- mean( kz.smp1 ) * mean( z.tbl1[,i][ h.lst1[[i]] ] )
    benefit.mnabove.seq2[i] <- mean( kz.smp2 ) * mean( z.tbl2[,i][ h.lst2[[i]] ] )
    benefit.mnabove.seq3[i] <- mean( kz.smp3 ) * mean( z.tbl3[,i][ h.lst3[[i]] ] )
    cost.mnabove.seq1   [i] <- mean( ka.smp1 ) * a.seq[i]
    cost.mnabove.seq2   [i] <- mean( ka.smp2 ) * a.seq[i]
    cost.mnabove.seq3   [i] <- mean( ka.smp3 ) * a.seq[i]
  }
  
  Rc.seq1  <- PRA.BDT1["R",] + cost.mnabove.seq1 - benefit.mnabove.seq1
  Rc.seq2  <- PRA.BDT2["R",] + cost.mnabove.seq2 - benefit.mnabove.seq2
  Rc.seq3  <- PRA.BDT3["R",] + cost.mnabove.seq3 - benefit.mnabove.seq3
  imaxRc.1 <- which( Rc.seq1 == min(Rc.seq1,na.rm=T) ) ; amaxRc.1 <- a.seq[imaxRc.1]
  imaxRc.2 <- which( Rc.seq2 == min(Rc.seq2,na.rm=T) ) ; amaxRc.2 <- a.seq[imaxRc.2]
  imaxRc.3 <- which( Rc.seq3 == min(Rc.seq3,na.rm=T) ) ; amaxRc.3 <- a.seq[imaxRc.3]
  
  par( mfrow=c(1,2), mar=c(4,2,3,0) )
  Rlim   <- range( 0, PRA.BDT1["R" ,], PRA.BDT2["R" ,], PRA.BDT3["R" ,], na.rm=T )
  Rclim  <- range( 0, Rc.seq1, Rc.seq2, Rc.seq3, na.rm=T )
  
  plot  ( a.seq, PRA.BDT1["R",], col="red"  , type="l",
          main=expression('R'), xlab="a", ylab="", xlim=c(0,3.5), ylim=Rlim )
  points( a.seq, PRA.BDT2["R",], col="red"  , type="l", lty=2 )
  points( a.seq, PRA.BDT3["R",], col="red"  , type="l", lty=3 )
  legend( "topright", legend=c("High unc","", "Low unc"), lty=3:1, cex=0.75, col="red" )
  
  plot  ( a.seq, Rc.seq1, col="red"  , type="l",
          main=expression( 'R'[c] ), xlab="a", ylab="", xlim=c(0,3.5), ylim=Rclim )
  points( a.seq, Rc.seq2, col="red"  , type="l", lty=2 )
  points( a.seq, Rc.seq3, col="red"  , type="l", lty=3 )
  abline( v=amaxRc.1    , col="black", lty=1 )
  abline( v=amaxRc.2    , col="black", lty=2 )
  abline( v=amaxRc.3    , col="black", lty=3 )
  text( amaxRc.1, min(Rc.seq1,na.rm=T), labels=paste("a_opt =",amaxRc.1 ), pos=4 )
  text( amaxRc.2, min(Rc.seq2,na.rm=T), labels=paste("a_opt =",amaxRc.2 ), pos=4 )
  text( amaxRc.3, min(Rc.seq3,na.rm=T), labels=paste("a_opt =",amaxRc.3 ), pos=4 )
  legend( "topright", legend=c("High unc","", "Low unc"), lty=3:1, cex=0.75, col="red" )
  
  ## Simplified Accounting for Cost and Benefit of the Action: $R_b$
  
  z.tbl1 <- z.tbl2 <- z.tbl3 <- NULL
  for(i in 1:np) {
    z.tbl1 <- rbind( z.tbl1, fz(a.seq,x.smp1[i],t.smp1[i]) + e.smp1[i] )
    z.tbl2 <- rbind( z.tbl2, fz(a.seq,x.smp2[i],t.smp2[i]) + e.smp2[i] )
    z.tbl3 <- rbind( z.tbl3, fz(a.seq,x.smp3[i],t.smp3[i]) + e.smp3[i] )
  }
  
  zmn.seq1 <- colMeans(z.tbl1)
  zmn.seq2 <- colMeans(z.tbl2)
  zmn.seq3 <- colMeans(z.tbl3)
  Rb.seq1  <- PRA.BDT1["R",] + mean(ka.smp1) * a.seq - mean(kz.smp1) * zmn.seq1
  Rb.seq2  <- PRA.BDT2["R",] + mean(ka.smp2) * a.seq - mean(kz.smp2) * zmn.seq2
  Rb.seq3  <- PRA.BDT3["R",] + mean(ka.smp3) * a.seq - mean(kz.smp3) * zmn.seq3
  imaxRb.1 <- which( Rb.seq1 == min(Rb.seq1,na.rm=T) ) ; amaxRb.1 <- a.seq[imaxRb.1]
  imaxRb.2 <- which( Rb.seq2 == min(Rb.seq2,na.rm=T) ) ; amaxRb.2 <- a.seq[imaxRb.2]
  imaxRb.3 <- which( Rb.seq3 == min(Rb.seq3,na.rm=T) ) ; amaxRb.3 <- a.seq[imaxRb.3]
  
  ## Only Correcting for Cost: $R_a$
  
  Ra.seq1  <- PRA.BDT1["R",] + a.seq * mean(ka.smp1)
  Ra.seq2  <- PRA.BDT2["R",] + a.seq * mean(ka.smp2)
  Ra.seq3  <- PRA.BDT3["R",] + a.seq * mean(ka.smp3)
  imaxRa.1 <- which( Ra.seq1 == min(Ra.seq1,na.rm=T) ) ; amaxRa.1 <- a.seq[imaxRa.1]
  imaxRa.2 <- which( Ra.seq2 == min(Ra.seq2,na.rm=T) ) ; amaxRa.2 <- a.seq[imaxRa.2]
  imaxRa.3 <- which( Ra.seq3 == min(Ra.seq3,na.rm=T) ) ; amaxRa.3 <- a.seq[imaxRa.3]
  
  par( mfrow=c(1,2), mar=c(4,2,3,0) )
  
  Rblim  <- range( 0, Rb.seq1, Rb.seq2, Rb.seq3, na.rm=T )
  
  plot  ( a.seq, Rb.seq1, col="red"  , type="l",
          main=expression('R'[b]), xlab="a", ylab="", xlim=c(0,3.5), ylim=Rblim )
  points( a.seq, Rb.seq2, col="red"  , type="l", lty=2 )
  points( a.seq, Rb.seq3, col="red"  , type="l", lty=3 )
  abline( v=amaxRb.1    , col="black", lty=1 )
  abline( v=amaxRb.2    , col="black", lty=2 )
  abline( v=amaxRb.3    , col="black", lty=3 )
  text( amaxRb.1, min(Rb.seq1,na.rm=T), labels=paste("a_opt =",amaxRb.1 ), pos=4 )
  text( amaxRb.2, min(Rb.seq2,na.rm=T), labels=paste("a_opt =",amaxRb.2 ), pos=4 )
  text( amaxRb.3, min(Rb.seq3,na.rm=T), labels=paste("a_opt =",amaxRb.3 ), pos=4 )
  legend( "topright", legend=c("High unc","","Low unc"), lty=3:1, cex=0.75, col="red" )
  
  Ralim  <- range( 0, Ra.seq1, Ra.seq2, Ra.seq3, na.rm=T )
  
  plot  ( a.seq, Ra.seq1, col="red"  , type="l",
          main=expression( 'R'[a] ), xlab="a", ylab="", xlim=c(0,3.5), ylim=Ralim )
  points( a.seq, Ra.seq2, col="red"  , type="l", lty=2 )
  points( a.seq, Ra.seq3, col="red"  , type="l", lty=3 )
  abline( v=amaxRa.1    , col="black", lty=1 )
  abline( v=amaxRa.2    , col="black", lty=2 )
  abline( v=amaxRa.3    , col="black", lty=3 )
  text( amaxRa.1, min(Ra.seq1,na.rm=T), labels=paste("a_opt =",amaxRa.1 ), pos=4 )
  text( amaxRa.2, min(Ra.seq2,na.rm=T), labels=paste("a_opt =",amaxRa.2 ), pos=4 )
  text( amaxRa.3, min(Ra.seq3,na.rm=T), labels=paste("a_opt =",amaxRa.3 ), pos=4 )
  legend( "topright", legend=c("High unc","","Low unc"), lty=3:1, cex=0.75, col="red" )
  