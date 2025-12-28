##############################



###   matrix



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





test_G <- function(m = 10,  p = 2, n = 100 ,P = 10){
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
    M_obs=array(0,c(m,m,n))
    
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
  
  
  U_true <- array(NA, c(n, P, P))
  for(i in 1:n){
    U_true[i,,] <- M_obs[,,i]
  }
  
  return(list(xx = X_obs, U_true = U_true))
}

#########################################################


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


simulation_matrix <- function(o_d = 50, out = 0.1,n = 100,q = 8){
  R_R <- c()
  NR_R <- c()
  for(e in 1:130){
    set.seed(12*e + 45)
    print(paste("+++",e))
    
    
    train_data <- data_G( p = 2, n = n ,P = q, out = out,o_d = o_d)
    test_data <- test_G( p = 2, n = 1 ,P = q)
    train_data$out_ind
    # BIC
    BIC_W <- function(Y, X, lam1, lam2){
      Y_hat <- array(NA,c(nrow(Y[,,1]),ncol(Y[1,,]),nrow(Y[1,,])))
      weight_i <- c()
      for(l in 1:nrow(Y[,,1])){
        model_bic <- RFRW(Y, X, X[l,], lam1, lam2)
        Y_hat[l,,] <- model_bic$U
        weight_i[l] <- model_bic$W[l]
      }
      
      dis_eta <- c()
      for(l in 1:nrow(Y[,,1])){
        dis_eta[l] <- DT(Y_hat[l,,], Y[l,,])
      }
      sigma_hat <- sum((weight_i*(dis_eta)^2))
      BIC <- length(dis_eta)*log((sigma_hat/(sum(weight_i)))) + (sum(weight_i < 1)) * log(length(dis_eta)+1)
      return(BIC = ifelse(sum(weight_i < 1) > length(dis_eta)*0.3, NA, BIC))
    }
    
    
    
    
    # Find the # Find the # Find the best parameter
    
    lam_set <- 400 * (seq(0.0000001, 1, length=20))^(0.8)
    lam_set <- as.matrix(expand.grid(lam_set, lam_set))
    L <- dim(lam_set)[1]
    test_mse_R <- c()
    test_mse_NR <- c()
    for(i in 1:nrow(test_data$xx)){
      print(i)
      res <- NULL
      for(l in 1:L){
        res[l] <- BIC_W(Y=train_data$Y, X=train_data$X, lam1=lam_set[l,1], lam2=lam_set[l,2])
      }
      rfit <- RFRW(train_data$Y, train_data$X, test_data$xx[i,], lam1=lam_set[which.min(res),1], lam2=lam_set[which.min(res),2])
      test_mse_R[i] <- DT(rfit$U, test_data$U_true[i,,])
      fit <- FR(Y = train_data$Y, X = train_data$X, xx = test_data$xx[i,])
      test_mse_NR[i] <- DT(fit$U, test_data$U_true[i,,])
      
    }
    
    R_R[e] <- mean(test_mse_R,na.rm=FALSE)
    R_R[e] <- ifelse(R_R[e] >=100, NA,R_R[e])
    NR_R[e] <- mean(test_mse_NR,na.rm=FALSE)
    NR_R[e] <- ifelse(R_R[e] >=100, NA,NR_R[e])
    
  }
  return(cbind(na.omit(NR_R),na.omit(R_R)))
}


n_can <- c(50,100)
o_d_can <- c(50,100)
out_can <- c(0,0.1,0.2)
q_can <- c(10)
data.var <- as.data.frame(expand.grid(n_can, o_d_can, out_can,q_can))
names(data.var) <- c("n_can", "o_d_can", "out_can","q_can")



pall_s_matrix <- function(qw){
  aa <- data.var[qw, "n_can"]
  bb <- data.var[qw, "o_d_can"]
  cc <- data.var[qw, "out_can"]
  dd <- data.var[qw, "q_can"]
  result_qw <- simulation_matrix(o_d = bb, out = cc,n = aa,q = dd)
  result_qw <- data.frame(result_qw)
  names(result_qw) <- c("NR_R","R_R")
  filename <- paste0("matrix","_n_", aa, "_bias_", bb, "_rate_", cc,"dim_",dd,".csv")
  return(result_qw)
  write.csv(result_qw, filename, row.names = FALSE)
}

library(parallel)


results <- mclapply(
  1:12,
  function(i) {
    tryCatch({
      message(paste("正在处理的模拟编号:", i))
      pall_s_matrix(qw = i)
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
  dd <- data.var[qw, "q_can"]
  data_save <- data.frame(results[[qw]])
  names(data_save) <- c("NR_R","R_R")
  filename <- paste0("matrix","_n_", aa, "_bias_", bb, "_rate_","_rate_", cc,"_p_",dd,".csv")
  write.csv(data_save, filename, row.names = FALSE)
}





