##############################



###   distrubution



###############################
rm(list=ls())



## robust method
RFRW <- function(Y, X, xx, lam1=10, lam2=10, max_itr=1000){
  # preparation 
  n <- dim(Y)[1]
  qq <- dim(Y)[2]
  mu_x <- mean(X)
  Sigma_x <- var(X)
  InvS <- solve(Sigma_x)
  
  # basic functions
  DT <- function(L_1, L_2){
    return( sqrt(sum((L_1 - L_2)^2)) + 10^(-10) )
  }
  
  # functions (weight for regression)
  Weight_x <- function(z, x){
    s2 <- as.vector( t(z-mu_x)%*%InvS%*%(x-mu_x) )
    return(1+s2)
  }
  
  # function (weight for outlier elimination)
  Weight_outlier <- function(r2, lam1, lam2) {
    if(r2<=lam1){
      return(1)
    }else if(r2>lam1 && r2<lam1+2*lam2){
      return(  1 - (r2-lam1)/(2*lam2) )
    }else{
      return(0)
    }
  }
  
  # function (update of matrix parameter)
  Update_U <- function(w_reg, w_out, Y){
    mat <- rep(0, ncol(Y))
    for(i in 1:n){
      mat <- mat + w_reg[i]*w_out[i]*Y[i,]
    }
    return( mat/sum(w_reg*w_out) )
  }
  
  # weight values for regression 
  w_reg <- rep(NA, n)
  for(i in 1:n){
    w_reg[i] <- Weight_x(xx, X[i])
  }
  
  # initial value 
  eta <- rep(0, n)    # distance-shift parameter
  U <- Y[which.min(abs(X - xx)),] + 0.00001
  
  
  # Iteration 
  for(itr in 1:max_itr){
    w_out <- c()
    for(i in 1:n){
      R2 <- w_reg[i]*(DT(Y[i,], U))^2
      w_out[i] <- Weight_outlier(R2, lam1, lam2)
    }
    U_new <- Update_U(w_reg, w_out, Y)
    dif_value <- (DT(U_new, U))^2
    U <- U_new
    if(dif_value<10^(-5)){ break() }
  }
  
  # Summary 
  Result <- list(U=U, itr=itr,W = w_out)
  return(Result)
}




## non-robust-method
FR <- function(Y, X, xx, max_itr=100){
  # preparation 
  n <- dim(Y)[1]
  qq <- dim(Y)[2]
  mu_x <- mean(X)
  Sigma_x <- var(X)
  InvS <- solve(Sigma_x)
  
  # basic functions 
  DT <- function(L_1, L_2){
    return( sqrt(sum((L_1 - L_2)^2)) + 10^(-10) )
  }
  
  # functions (weight for regression)
  Weight_x <- function(z, x){
    s2 <- as.vector( t(z-mu_x)%*%InvS%*%(x-mu_x) )
    return(1+s2)
  }
  
  # function (update of matrix parameter)
  Update_U <- function(w_reg, Y){
    mat <- rep(0, ncol(Y))
    for(i in 1:n){
      mat <- mat + w_reg[i]*Y[i,]
    }
    return( mat/sum(w_reg) )
  }
  
  # weight values for regression 
  w_reg <- rep(NA, n)
  for(i in 1:n){
    w_reg[i] <- Weight_x(xx, X[i])
  }
  
  # initial value 
  eta <- rep(0, n)    # distance-shift parameter
  U <- Y[which.min(abs(X - xx)),] + 0.00001
  
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

