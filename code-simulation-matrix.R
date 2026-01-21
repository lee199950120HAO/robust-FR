##############################



###   matrix



###############################



rm(list=ls())



tr <- function(x){ sum(diag(x)) }
DT <- function(L_1, L_2){
  mat <- t(L_1-L_2)%*%(L_1-L_2)
  return( sqrt(tr(mat))+10^(-10) )
}
## non-robust method
FR <- function(Y, X, xx, max_itr=100){
  # preparation
  n <- dim(Y)[1]
  qq <- dim(Y)[2]
  mu_x <- apply(X, 2, mean)
  Sigma_x <- var(X)
  InvS <- solve(Sigma_x)
  
  # basic functions
  tr <- function(x){ sum(diag(x)) }
  DT <- function(L_1, L_2){
    mat <- t(L_1-L_2)%*%(L_1-L_2)
    return( sqrt(tr(mat))+10^(-10) )
  }
  
  # functions (weight for regression)
  Weight_x <- function(z, x){
    s2 <- as.vector( t(z-mu_x)%*%InvS%*%(x-mu_x) )
    return(1+s2)
  }
  
  # function (update of matrix parameter)
  Update_U <- function(w_reg, Y){
    mat <- matrix(0, qq, qq)
    for(i in 1:n){
      mat <- mat + w_reg[i]*Y[i,,]
    }
    return( mat/sum(w_reg) )
  }
  
  # weight values for regression
  w_reg <- rep(NA, n)
  for(i in 1:n){
    w_reg[i] <- Weight_x(xx, X[i,])
  }
  
  # initial value
  sel <- which.min( colMeans((t(X)-xx)^2) )
  U <- Y[sel,,]   # initial value of the target matrix
  
  # Iteration
  for(itr in 1:max_itr){
    U_new <- Update_U(w_reg, Y)
    dif_value <- DT(U_new, U)
    U <- U_new
    if(dif_value<10^(-5)){ break() }
  }
  
  # Summary
  Result <- list(U=U, itr=itr)
  return(Result)
}



RFRW <- function(Y, X, xx, lam1=10, lam2=10, max_itr=100){
  # preparation
  n <- dim(Y)[1]
  qq <- dim(Y)[2]
  mu_x <- apply(X, 2, mean)
  Sigma_x <- var(X)
  InvS <- solve(Sigma_x)
  
  # basic functions
  tr <- function(x){ sum(diag(x)) }
  DT <- function(L_1, L_2){
    mat <- t(L_1-L_2)%*%(L_1-L_2)
    return( sqrt(tr(mat))+10^(-10) )
  }
  
  # functions (weight for regression)
  Weight_x <- function(z, x){
    s2 <- as.vector( t(z-mu_x)%*%InvS%*%(x-mu_x) )
    return(1+s2)
  }
  
  # function (weight for 0.1lier elimination)
  Weight_0.1lier <- function(r2, lam1, lam2) {
    if(r2<=lam1){
      return(1)
    }else if(r2>lam1 && r2<lam1+2*lam2){
      return(  1 - (r2-lam1)/(2*lam2) )
    }else{
      return(0)
    }
  }
  
  # function (update of matrix parameter)
  Update_U <- function(w_reg, w_0.1, Y){
    mat <- matrix(0, qq, qq)
    for(i in 1:n){
      mat <- mat + w_reg[i]*w_0.1[i]*Y[i,,]
    }
    return( mat/sum(w_reg*w_0.1) )
  }
  
  # weight values for regression
  w_reg <- rep(NA, n)
  for(i in 1:n){
    w_reg[i] <- Weight_x(xx, X[i,])
  }
  
  # initial value
  sel <- which.min( colMeans((t(X)-xx)^2) )
  U <- Y[sel,,]   # initial value of the target matrix
  
  # Iteration
  for(itr in 1:max_itr){
    w_0.1 <- c()
    for(i in 1:n){
      R2 <- w_reg[i]*DT(Y[i,,], U)
      w_0.1[i] <- Weight_0.1lier(R2, lam1, lam2)
    }
    U_new <- Update_U(w_reg, w_0.1, Y)
    dif_value <- DT(U_new, U)
    U <- U_new
    if(dif_value<10^(-5)){ break() }
  }
  
  # Summary
  Result <- list(U=U, itr=itr, W=w_0.1)
  return(Result)
}




