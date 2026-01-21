##############################



###   matrix response type-2



###############################



rm(list=ls())

library("frechet")
library("shapes")
library("stringr")
library("pbapply")
library("parallel")
library("MASS")
library("matrixcalc")






data_G <- function(p = 2, n = 100 ,P = 10, out = 0.1,o_d = 50){
  run <- 100  #number of main loops
  ISE <- matrix(0,run,4)
  matrix.exp <- function(A){
    eig <- eigen(A)
    EA=eig$vectors%*%diag(exp(eig$values))%*%t(eig$vectors)
    return((EA+t(EA))/2)
  }
  ##Model settings################################################
  setting <- 2
  
  method <- "Log_cholesky"; metric <- "log_cholesky"; alpha <- 1 #type of metric
  if (setting ==1) {
    if (p==2) {beta <- c(0.75,0.25)}
    if (p==5) {beta <- c(0.1,0.2,0.3,0.4,0)}
    if (p==10) {beta <- c(0.1,0.2,0.3,0.4,rep(0,6))}
    if (p==20) {beta <- c(0.1,0.2,0.3,0.4,rep(0,12),0.1,0.2,0.3,0.4)/2}
  }
  if (setting ==2) {
    if (p==2) {beta_1 <- c(0.75,0.25);beta_2 <- c(0.25,0.75)}
    if (p==5) {beta_1 <- c(0.1,0.2,0.3,0.4,0);beta_2 <- c(0,0.1,0.2,0.3,0.4)}
    if (p==10) {beta_1 <- c(0.1,0.2,0.3,0.4,rep(0,6));beta_2 <- c(rep(0,6),0.1,0.2,0.3,0.4)}
    if (p==20) {beta_1 <- c(0.1,0.2,0.3,0.4,rep(0,16));beta_2 <- c(rep(0,16),0.1,0.2,0.3,0.4)}
  }
  
  for (bb in 1:run) {
    data <- runif(n*p, 0, 1)
    X_obs <- matrix(data, n, p) #trainning data
    M_obs=array(0,c(P,P,n))
    
    for(j in 1:n){
      pho <- c()
      for(i in 1:(P-1)){
        pho[i] <- round(runif(1),1)*cos(as.numeric(t(beta_1)%*%X_obs[j,])*4*pi)
      }
      
      D <- matrix(0,P,P)
      for(i in 1:P){
        D[,i] <- c(1:P-(i-1))
        
      }
      D <- D - 1
      DD <- matrix(0,P,P)
      for(i in 1:(P-1)){
        for(k in 1:nrow(which(D == i,arr.ind = TRUE))){
          DD[which(D == i,arr.ind = TRUE)[k,][1],which(D == i,arr.ind = TRUE)[k,][2]] <- pho[i]
        }
      }
      D <- t(DD) + DD + diag(P)
      
      Z <- matrix(0,P,P)
      for (i in 1:P) {
        Z[i, i] <- rnorm(1, 0, 1)
      }
      for (i in 1:(P - 1)) {
        for (k in (i + 1):P) {
          val <- rnorm(1, 0, 1 / sqrt(2))
          Z[i, k] <- val
          Z[k, i] <- val
        }
      }
      logY <- 0.2*Z+D
      
      M_obs[,,j] <-matrix.exp(logY)
    }
    
  }

  out_ind <- rep(0, n)
  out_ind[sample(1:n, round(n*out))] <- 1
  for(i in 1:n){
    if(out_ind[i]==1){
      M_obs[,,i] <- M_obs[,,i] + o_d*matrix(rbinom(P^2, 1, 1), P, P)
    }
  }


  U_true <- array(NA, c(n, P, P))
  for(i in 1:n){
    U_true[i,,] <- M_obs[,,i]
  }
  return(list(X = X_obs, Y = U_true,out_ind = out_ind,M_obs = M_obs))
}

