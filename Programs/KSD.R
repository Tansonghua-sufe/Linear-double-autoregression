# Note: This file contains basic functions of KSD test.
#       The method is based on Luo, D., Zhu, K., Gong, H., & Li, D. (2021). Testing error distribution by kernelized Stein discrepancy in multivariate time series models. Journal of Business & Economic Statistics, 1-15.
#       The code is a modified version of R package KSD: Goodness-of-Fit Tests using Kernelized Stein Discrepancy
#       see https://cran.r-project.org/web/packages/KSD/

score_f_R <- function(x,p="norm",df.t=5){
  # Score funtion
  return(switch (p,
                 "norm" = -x,
                 "t"=-(df.t+1)*x/(df.t-2+x^2),
                 "laplace" = -sqrt(2)*ifelse(x>0,1,-1)
  ))
}

find_median_distance <- function(Z){
  # find the median of residual distance
  if(is.data.frame(Z)){
    Z = data.matrix(Z)
  }else{
    Z = as.array(Z)
  }
  size1 <- dim(Z)[1]
  size2 <- dim(Z)[2]

  Zmed <- Z

  Zmedsq <- Zmed * Zmed;
  if(is.na(dim(Z)[2]))
    G <- Zmedsq
  else
    G <- rowSums(Zmedsq)

  # Create row/col repeated matrices
  Q <- rep.col(G,size1)
  # Q <- matrix(rep(G,size1),nrow = length(G),ncol = size1)
  R <- rep.row(t(G),size1)
  # R <- t(Q)

  dists <- Q + R - 2 * Zmed %*% t(Zmed)
  dists[lower.tri(dists, diag = TRUE)] = 0
  dists <- array(dists,dim=c(size1^2,1))
  median_dist <- median(dists[dists > 0 ])

  return(median_dist)
}

rep.col<-function(x,n){
  # auxiliary function
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

rep.row<-function(x,n){
  # auxiliary function
  matrix(rep(x,each=n),nrow=n)
}

repmat <- function(X,a=1,b=1){
  # auxiliary function
  rows <- dim(X)[1]
  cols <- dim(X)[2]
  if(is.null(cols))
    cols <- 1
  rowRep <- matrix(rep(t(X),a),ncol = cols, byrow = TRUE)
  newX <- matrix(rep(rowRep, b),ncol=cols*b)
  return(newX)
}

ksd_R <- function(x,p,sigma){

  n <- length(x)
  h <- sigma

  Sqx = score_f_R(x,p)

  XY <- x %*% t(x)
  x2 <- x^2
  sumx <- x * Sqx

  X2e <- repmat(x2,1,n)
  H <- (X2e + t(X2e) - 2 * XY)
  Kxy <- exp(-H/(2 * h^2))
  sqxdy <- -(Sqx %*% t(x) - repmat(sumx,1,n))/h^2
  dxsqy <- t(sqxdy)
  dxdy <- (-H/h^4 + 1/h^2)
  M = (Sqx %*% t(Sqx) + sqxdy + dxsqy + dxdy) * Kxy
  M2 <- M - diag(diag(M))
  return(sum(M2)/(n * (n - 1)))
}



bootstrap_KSD <- function(lambda_est,distribution = "laplace",p,N,B,est_method = "G-QMLE"){

  # ==========================================================
  # Bootstrap function of KSD
  # lambda_est: estimated parameter
  # distribution: the distribution we want to test
  # p: order
  # N: sample size
  # B: bootstrap sample size
  # est_method: the method of estimation ("G-QMLE" or "E-QMLE")
  # ==========================================================
  burning <- 200

  record <- array(dim = c(B))
  for(b in 1:B){
    # generating new innovations based on different distribution assumption
    eta_b <- switch (distribution,
                     "laplace" = rlaplace(burning+N),
                     "norm" = rnorm(burning+N)
    )
    # standardization based on different method of estimation
    eta_b <- switch (est_method,
                     "G-QMLE" = eta_b/sd(eta_b),
                     "E-QMLE" = eta_b/mean(abs(eta_b))
    )
    y_b <- DGP(epsilon = eta_b,lambda_est[1:p],lambda_est[(p+1):(2*p+1)]) # Generating simulated data by given parameter
    Y_b <- y_b[(burning+1):(burning+N)] # only the last N observations

    inf_bound <- c(rep(-1,p),0.01,rep(1e-8,p))
    sup_bound <- c(rep(1,p),20,rep(20,p))

    fit_b<- switch (est_method, # estimation
                    "G-QMLE" = optim(par = lambda_est ,fn=likehood_LDAR_GQMLE,lower= inf_bound,upper = sup_bound,method = "L-BFGS-B",Y=Y_b),
                    "E-QMLE" = optim(par = lambda_est ,fn=likehood_LDAR_EQMLE,lower= inf_bound,upper = sup_bound,method = "L-BFGS-B",Y=Y_b)
    )

    lambda_b <- fit_b$par
    res_b <- resdual(Y_b,lambda_b) # residuals

    res_s_b <- (res_b-mean(res_b))/sd(res_b) # standardization

    sigma_b <- find_median_distance(Z = res_s_b)
    test_b <- ksd_R(x = res_s_b,p = distribution,sigma = sigma_b) # test statistic
    record[b] <- test_b
  }
  return(record)
}
