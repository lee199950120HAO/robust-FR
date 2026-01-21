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
  
  Weight_outlier <- function(r2, lam1, lam2) {
    if(r2<=lam1){
      return(1)
    }else if(r2>lam1 && r2<lam1+2*lam2){
      return(  1 - (r2-lam1)/(2*lam2) )
    }else{
      return(0)
    }
  }
  

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
  n <- dim(Y)[1]
  qq <- dim(Y)[2]
  mu_x <- mean(X)
  Sigma_x <- var(X)
  InvS <- solve(Sigma_x)
  
  DT <- function(L_1, L_2){
    return( sqrt(sum((L_1 - L_2)^2)) + 10^(-10) )
  }
  
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
  
  w_reg <- rep(NA, n)
  for(i in 1:n){
    w_reg[i] <- Weight_x(xx, X[i])
  }
  
  # initial value 
  eta <- rep(0, n)    
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







simulation_distrubution <- function(n = 100, o_d = 100, out = 0.1, X_test = 0.5){
  R_R <- c()
  NR_R <- c()
  for(e in 1:300){
    set.seed(12*e + 45)
    print(paste("+++",e))
    mu0 <- 0
    sigma0 <- 3
    beta <- 3
    gamma <- 0.5
    v1 <- 0.25
    v2 <- 1
    # covariate
    X <- runif(n, 0, 1)
    
    # true mu and sigma
    mu <- rnorm(n, mean = mu0 + beta * X, sd = sqrt(v1))
    sigma <- rgamma(n, shape = (sigma0 + gamma * X)^2 / v2, rate = v2 / (sigma0 + gamma * X))
    
    
    # generation of response
    z <- seq(0.1, 0.9, by=0.01)
    Y <- matrix(0, ncol=length(z), nrow=n)
    
    out_ind <- rep(1, nrow(Y))
    out_ind[sample(1:n, round(n*out))] <- 0
    for(i in 1:n){
      for(j in 1:length(z)){
        Y[i,j] <- mu[i] + sigma[i]*qnorm(z[j],0,1)  
      }
      if(out_ind[i]==0){
        Y[i,] <- Y[i,] + o_d
      }
    }
    
    
    matplot(t(Y), type="l", col="black", lty=1)
    xx <- X_test
    mu_xx <- rnorm(1, mean = mu0 + beta * xx, sd = sqrt(v1))
    sigma_xx <- rgamma(1, shape = (sigma0 + gamma * xx)^2 / v2, rate = v2 / (sigma0 + gamma * X))
    
    Y_test <- c()
    for(i in 1:length(z)){
      Y_test[i] <- mu_xx + sigma_xx*qnorm(z[i])
    }
    points(Y_test, col="red", type="l", lty=1, lwd=2)
    
    
    ##distance
    distance <- function(x,y){
      return(sqrt(sum((x - y)^2)))
    }
    
    xx <- X_test
    xx_est <- Y_test
    # robust method 
    
    lam_set <- 82000 * (seq(0.0000001, 1, length=20))^(0.8)
    lam_set <- as.matrix(expand.grid(lam_set, lam_set))
    L <- dim(lam_set)[1]
    
   # BIC
  BIC_W <- function(Y, X, lam1, lam2){
    Y_hat <- matrix(data = NA, nrow = nrow(Y),ncol = ncol(Y))
    weight_i <- c()
    for(l in 1:nrow(Y)){
      model_bic <- RFRW(Y, X, X[l], lam1 = lam1, lam2 = lam2)
      Y_hat[l,] <- model_bic$U
      weight_i[l] <- model_bic$W[l]
    }
    
    dis_eta <- c()
    for(l in 1:nrow(Y)){
      dis_eta[l] <- distance(Y_hat[l,], Y[l,])
    }
    
    
    sigma_hat <- sum((weight_i*(dis_eta)^2))
    BIC <- length(dis_eta)*log((sigma_hat/(sum(weight_i)))) + (sum(weight_i < 1)) * log(length(dis_eta)+1)
    return(BIC = ifelse(sum(weight_i < 1) > length(dis_eta)*0.3, NA, BIC))
  }
    
    # Find the best parameter
    res <- NULL
    for(l in 1:L){
      res[l] <- BIC_W(Y=Y, X=X, lam1=lam_set[l,1], lam2=lam_set[l,2])
      print(l)
    }
    
    
    rfit <- RFRW(Y, X, xx, lam1=lam_set[which.min(res),1], lam2=lam_set[which.min(res),2])
    
    
    
    # non-robust method
    fit <- FR(Y, X, xx)
    R_R[e] <- distance(rfit$U, xx_est)
    NR_R[e] <- distance(fit$U, xx_est)
    if(length(R_R) == 100) break
  }
  return(cbind(na.omit(NR_R),na.omit(R_R)))
}



n_can <- c(50,100)
o_d_can <- c(50,100)
out_can <- c(0,0.1,0.2)

data.var <- as.data.frame(expand.grid(n_can, o_d_can, out_can))
names(data.var) <- c("n_can", "o_d_can", "out_can")

pall_s_distrubution <- function(qw){
  aa <- data.var[qw, "n_can"]
  bb <- data.var[qw, "o_d_can"]
  cc <- data.var[qw, "out_can"]
  result_qw <- simulation_distrubution(n = aa, o_d = bb, out = cc, X_test = 0.5)
  names(result_qw) <- c("NR_R","R_R") 
  filename <- paste0("distrubution","_n_", aa, "_bias_", bb, "_rate_", cc,".csv")
  return(result_qw)
  write.csv(result_qw, filename, row.names = FALSE)
}
library(parallel)

results <- mclapply(
  1:12,
  function(i) {
    tryCatch({
      pall_s_distrubution(qw = i)
    }, error = function(e) {
      cat("simulation: ", i, "failed: ",conditionMessage(e) )
    }
    )
    
  },
  mc.cores=9
)


for(qw in 1:12){
  aa <- data.var[qw, "n_can"]
  bb <- data.var[qw, "o_d_can"]
  cc <- data.var[qw, "out_can"]
  data_save <- data.frame(results[[qw]])
  names(data_save) <- c("NR_R","R_R") 
  filename <- paste0("distrubution","_n_", aa, "_bias_", bb, "_rate_", cc,".csv")
  write.csv(data_save, filename, row.names = FALSE)
}


