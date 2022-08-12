estimate_EQMLE<- function(Y,p){
  ## Parameter estimation of LDAR-EQMLE
  # ==========Input====================
  # Y: Observation
  # p: Order
  # ===================================
  init_value <- c(rep(0.5,p),1,rep(0.5,p))
  inf_bound <- c(rep(-1,p),0.01,rep(1e-8,p))
  sup_bound <- c(rep(20,p),20,rep(20,p))
  fit_EQMLE<- optim(par = init_value ,fn=likehood_LDAR_EQMLE,lower= inf_bound,upper = sup_bound,method = "L-BFGS-B",Y=Y)
  lambda_EQMLE <- fit_EQMLE$par
  names(lambda_EQMLE) <- c(paste("alpha",1:p,sep = ""),"omega",paste("beta",1:p,sep = ""))
  return(lambda_EQMLE)
}


asd_EQMLE <- function(Y,p,lambda=NULL,cov=F){
  ## ASD or covariance of LDAR-EQMLE
  # ==========Input=================================================================
  # Y: Observation
  # p: Order
  # lambda: Parameter estimation (if NULL, use `estimate_EQMLE` to obtain the value)
  # cov: If cov==F, ASD is return; If cov==T, covariance is return
  # ================================================================================
  N <- length(Y)
  if(is.null(lambda)){
    lambda <- estimate_EQMLE(Y,p)
  }

  res <- resdual(Y,lambda)
  E_1 <- mean(res)
  E_2 <- mean(res^2)
  n_0 <- length(res)-p # effective sample size
  IQR <- quantile(res, 0.75) - quantile(res, 0.25)
  a <- min(sd(res),IQR/1.34)
  h <- 0.9*n_0^(-1/5)*a
  d <- stats::density(res, bw = h, kernel = "gaussian")
  x_negative <- sum(d$x<0);x_positive <- x_negative+1
  f_0 <- (d$y[x_negative]+d$y[x_positive])/2
  covariance <- Xi_EQMLE(lambda,Y,E_1,E_2,f_0)
  if(cov==F){
    ASD <- sqrt(diag(covariance)/N)
    names(ASD) <- c(paste("alpha",1:p,sep = ""),"omega",paste("beta",1:p,sep = ""))
    return(ASD)
  }else{
    return(covariance)
  }
}

estimate_GQMLE<- function(Y,p){
  ## Parameter estimation of LDAR-GQMLE
  # ==========Input====================
  # Y: Observation
  # p: Order
  # ===================================
  init_value <- c(rep(1,p),1,rep(1,p))
  inf_bound <- c(rep(-1,p),0.01,rep(1e-8,p))
  sup_bound <- c(rep(20,p),20,rep(20,p))
  fit_GQMLE<- optim(par = init_value ,fn=likehood_LDAR_GQMLE,
                   lower= inf_bound,upper = sup_bound,method = "L-BFGS-B",Y=Y)
  lambda_GQMLE <- fit_GQMLE$par
  names(lambda_GQMLE) <- c(paste("alpha",1:p,sep = ""),"omega",paste("beta",1:p,sep = ""))
  return(lambda_GQMLE)
}

asd_GQMLE <- function(Y,p,lambda=NULL,cov=F){
  ## ASD or covariance of LDAR-GQMLE
  # ==========Input=================================================================
  # Y: Observation
  # p: Order
  # lambda: Parameter estimation (if NULL, use `estimate_GQMLE` to obtain the value)
  # cov: If cov==F, ASD is return; If cov==T, covariance is return
  # ================================================================================
  N <- length(Y)
  if(is.null(lambda)){
    lambda <- estimate_GQMLE(Y,p)
  }
  res <- resdual(Y,lambda_GQMLE)
  E_3 <- mean(res^3)
  E_4 <- mean(res^4)
  covariance <- Xi_GQMLE(lambda_GQMLE,Y,E_3 = E_3,E_4 = E_4)
  if(cov==F){
    ASD <- sqrt(diag(covariance)/N)
    names(ASD) <- c(paste("alpha",1:p,sep = ""),"omega",paste("beta",1:p,sep = ""))
    return(ASD)
  }else{
    return(covariance)
  }
}

modelsection_LDAR_EQMLE <- function(Y,P_max){
  ## Model section (Based on EQMLE)
  # ==========Input======================================
  # Y: Observation
  # P_max: The max lag, predetermined positive integer.
  # =====================================================
  n=length(Y)
  recommended_order<-vector(length = 1) # the recommended order
  BIC_min<-1e+20
  for (p in 1:P_max){
    lambda_EQMLE <- estimate_GQMLE(Y = Y,p = p)
    likelihood_BIC_EQMLE <- likehood_EQMLE_p_max(lambda_EQMLE,Y,P_max) # likelihood
    BIC <- 2*likelihood_BIC_EQMLE+(2*p+1)*log(n-P_max) # BIC
    if (BIC<BIC_min){
      recommended_order[1] <- p
      BIC_min <- BIC
    }
  }
  return(recommended_order)
}

predict_EQMLE<-function(Y_train,Y_test,p,q,method = "rolling"){
  ## One-step rolling forecasting based on LDAR-EQMLE
  # ==========Input======================================
  # Y_train: Train set
  # Y_test: Test set
  # p: Order
  # q: Interested quantile levels
  # method: rolling forecasting with rolling start point ("rolling")
  #         or fixed start point ("fixed")
  # =====================================================
  n_train<-length(Y_train)
  n_test <- length(Y_test)
  Y_hat_mean <- array(dim = c(length(Y_test)))
  Y_hat_q <- array(dim = c(length(Y_test),length(q)))
  for (t in 1:n_test) {
    if(method == "rolling"){Y <- c(Y_train[t:n_train],Y_test[0:(t-1)])} #fixed window
    if(method == "fixed"){Y <- c(Y_train[1:n_train],Y_test[0:(t-1)])} #fixed start
    n <- length(Y)
    Y_1 <- rev(Y[n:(n-p+1)])
    Y_2 <- abs(Y_1)
    lambda <- estimate_EQMLE(Y,p)
    res <- resdual(Y,lambda)
    lamba_1<-lambda[1:p]
    omega <- lambda[p+1]
    lambda_2 <- lambda[(p+2):(2*p+1)]
    Y_hat_mean[t] <- Y_1%*%lamba_1
    Y_hat_q[t,] <- Y_hat_mean[t]+quantile(x=res,probs = q)*(omega+Y_2%*%lambda_2)
  }
  return(list('Y_hat_mean'=Y_hat_mean,"Y_hat_q"=Y_hat_q))
}


predict_GQMLE<-function(Y_train,Y_test,q,p,method = "rolling"){
  ## One-step rolling forecasting based on LDAR-GQMLE
  # ==========Input======================================
  # Y_train: Train set
  # Y_test: Test set
  # p: Order
  # q: Interested quantile levels
  # method: rolling forecasting with rolling start point ("rolling")
  #         or fixed start point ("fixed")
  # =====================================================
  n_train <- length(Y_train)
  n_test <- length(Y_test)
  Y_hat_mean <- array(dim = c(length(Y_test)))
  Y_hat_q <- array(dim = c(length(Y_test),length(q)))

  for (t in 1:n_test) {
    if(method == "rolling"){Y <- c(Y_train[t:n_train],Y_test[0:(t-1)])}
    if(method == "fixed"){Y <- c(Y_train[1:n_train],Y_test[0:(t-1)])}
    n <- length(Y)
    Y_1 <- rev(Y[n:(n-p+1)])
    Y_2 <- abs(Y_1)
    fitted_value <- estimate_GQMLE(Y,p)
    res <- resdual(Y,fitted_value)
    lambda_1 <- fitted_value[1:p]
    omega <- fitted_value[p+1]
    lambda_2 <- fitted_value[(p+2):(2*p+1)]

    Y_hat_mean[t] <- Y_1%*%lambda_1
    Y_hat_q[t,] <- Y_hat_mean[t]+quantile(x = res,probs = q)*(omega+Y_2%*%lambda_2)
  }
  return(list("Y_hat_mean"=Y_hat_mean,"Y_hat_q"=Y_hat_q))
}

predict_DWQRE<-function(Y_train,Y_test,q,p,method = "rolling"){
  # ==========Input======================================
  # Y_train: Train set
  # Y_test: Test set
  # p: Order
  # q: Interested quantile levels
  # method: rolling forecasting with rolling start point ("rolling")
  #         or fixed start point ("fixed")
  # =====================================================
  n_train <- length(Y_train)
  n_test <- length(Y_test)
  Y_hat_mean <- array(dim = c(length(Y_test)))
  Y_hat_q <- array(dim = c(length(Y_test),length(q)))

  for (t in 1:n_test) {
    if(method == "rolling"){Y <- c(Y_train[t:n_train],Y_test[0:(t-1)])}
    if(method == "fixed"){Y <- c(Y_train[1:n_train],Y_test[0:(t-1)])}
    n <- length(Y)
    Y_1 <- rev(Y[n:(n-p+1)])
    Y_2 <- abs(Y_1)
    ldar=ldarest(n,p,lagm,K,KK,ell,Y)
    residual=ldar$residual
    lambda_1 <-ldar$Phihat3
    lambda_2 <- ldar$Betahat3
    Y_hat_mean[t] <- Y_1%*%lambda_1
    Y_hat_q[t,] <- Y_hat_mean[t]+quantile(x = residual,probs = q)*(1+Y_2%*%lambda_2)
  }
  return(list("Y_hat_mean"=Y_hat_mean,"Y_hat_q"=Y_hat_q))
}



VaRbacktestB <- function(Hit,VaR,tau,p) {
  ## VaR backtest function
  # ============================================================
  # Hit is a sequence of Hit_t=I(y_t<-VaR_t)
  # VaR is a sequence of Value-at-Risk
  # tau is the quantile level of -VaR=Q_tau(y_t|F_{t-1})
  # p is the dimension of lagged hits in DQ tests
  # ============================================================
  n <- length(Hit) # sample size
  ### Uncond coverage test
  n1 <- sum(Hit); tauhat <- n1/n
  likUC_null <- tau^n1*(1-tau)^(n-n1); likUC_alt <- tauhat^n1*(1-tauhat)^(n-n1)
  LR_UC <- -2*log(likUC_null/likUC_alt)
  P.UC <- pchisq(LR_UC,df=1,lower.tail=FALSE) # p-value
  ### Independence test
  ContTable <- table(Hit[-n],Hit[-1])
  n00 <- ContTable[1,1]; n01 <- ContTable[1,2]; n10 <- ContTable[2,1]; n11 <- ContTable[2,2]
  tau_null <- (n01+n11)/(n00+n01+n10+n11)
  tau0_alt <- n01/(n00+n01); tau1_alt <- n11/(n10+n11)
  likInd_null <- tau_null^(n01+n11)*(1-tau_null)^(n00+n10)
  likInd_alt <- tau0_alt^n01*(1-tau0_alt)^n00*tau1_alt^n11*(1-tau1_alt)^n10
  LR_Ind <- -2*log(likInd_null/likInd_alt)
  P.Ind <- pchisq(LR_Ind,df=1,lower.tail=FALSE) # p-value
  ### Cond coverage test
  LR_CC <- LR_UC+LR_Ind
  P.CC <- pchisq(LR_CC,df=2,lower.tail=FALSE)
  ### DQ test 1: hits only
  X <- cbind(1,embed(Hit[-n],p)) # (n-p)*(p+1) matrix
  if (det(t(X)%*%X)<=10^{-5}){P.DQ1 <- NaN
  } else {
    DQ1 <- t(Hit[(p+1):n]-tau)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%(Hit[(p+1):n]-tau)/tau/(1-tau)
    P.DQ1 <- pchisq(DQ1,df=p+1,lower.tail=FALSE) # p-value
  }
  ### DQ test 2: hits & lagged VaR
  X <- cbind(1,embed(Hit[-n],p),VaR[(p+1):n]) # (n-p)*(p+2) matrix
  if (det(t(X)%*%X)<=10^{-5}){P.DQ2=NaN
  } else {
    DQ2 <- t(Hit[(p+1):n]-tau)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%(Hit[(p+1):n]-tau)/tau/(1-tau)
    P.DQ2 <- pchisq(DQ2,df=p+2,lower.tail=FALSE) # p-value
  }
  #   ### Absolute deviation of ECR: abs(ECR%-tau%)*100
  ECR <- mean(Hit)
  rbind.data.frame(
    ECR=ECR,UC=P.UC,Ind=P.Ind,CC=P.CC,DQ_hits=P.DQ1,DQ_comb=P.DQ2)  #,deparse.level=0
}
