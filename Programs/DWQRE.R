#########################################################################
#-----------------------------------------------------------------------#
# Created on:    2017/12/13                                             #
# Updated on:    2017/12/28                                             #
# EMAIL:         zhunanapig@126.com                                     #
#                SUFE/HKU, ZHU Qianqian. (2017)                         #
#########################################################################
library(quantreg)
library(MASS)  #mvrnorm
library(Matrix)  #nearPD

##########################################################################
# The following functions are used to calculate optimal weights and ASDs #
##########################################################################

density=function(btau,resid,N,h,K)
{
  #btau: K*1 vector of quantiles; resid: (n-lagm)*1 vector of residuals;
  #N: sample size; h: bandwidth; K: number of quantiles
  dens=NULL
  for (i in 1:K)
  {
    dens[i]=sum(dnorm((btau[i]-resid)/h))/(N*h)
  }
  return(dens)
}
gamma=function(K,Tau,density)
{
  #density:K*1 vector of density for Tau; Tau: K*1 vector of quantiles
  gamma=matrix(0,K,K)
  for (i in 1:K)
  {
    for (j in 1:K)
    {
      gamma[i,j]=(min(Tau[i],Tau[j])-Tau[i]*Tau[j])/(density[i]*density[j])
    }
  }
  gamma
}
asd=function(p,N,K,target,X,XXt,beta,phi,btau,Tau,wt)
{ ###obtain asd of beta and phi and optimal weight
  #pi_{k,opt}  #2p*2p matrix
  #V=(V_{kk'})  #2pK*2pK matrix
  ##########Gamma######
  resid=(target-X%*%phi)/(1+abs(X)%*%beta)  #residual
  h=0.9*min(as.numeric(quantile(resid,0.75)-quantile(resid,0.25))/1.34,sd(resid))*N^(-1/5)  #bandwidth
  dens=density(btau,resid,N,h,K)  #K*1 vector
  gamma=gamma(K,Tau,dens) #K*K matrix
  ##########Sigma1#########
  Sigma1=array(0,dim=c(2*p,2*p,K))
  for (i in 1:K)
  {
    Sigma1[,,i]=rbind(cbind(btau[i]*diag(p),matrix(0,p,p)),cbind(matrix(0,p,p),diag(p)))
  }
  ##########Omega0&Omega1&Omega2#########
  B=matrix(0,2*p+1,2*p+1)
  C=matrix(0,2*p+1,2*p+1)
  for (i in 1:(N-lagm))
  {
    B=B+wt[i]/as.numeric(1+abs(X[i,])%*%beta)*t(matrix(XXt[i,],nrow=1))%*%matrix(XXt[i,],nrow=1)
    C=C+wt[i]^2*t(matrix(XXt[i,],nrow=1))%*%matrix(XXt[i,],nrow=1)
  }
  Omega0=B/N
  Omega1=solve(Omega0)%*%(C/N)%*%solve(Omega0)
  D=rbind(cbind(-beta,diag(p),matrix(0,p,p)),cbind(matrix(0,p,p+1),diag(p)))
  Omega2=D%*%Omega1%*%t(D)
  ##########V#########
  V=matrix(0,2*p*K,2*p*K)
  for(i in 1:K)
  {
    for (j in 1:K)
    {
      V[((i-1)*2*p+1):(i*2*p),((j-1)*2*p+1):(j*2*p)]=gamma[i,j]*Omega2
    }
  }
  ##########Vs&Vss#########
  Vs=array(0,dim=c(2*p,2*p,K))
  sigma=array(0,dim=c(2*p,2*p,K))
  for (k in 1:K)
  {
    for (i in 1:K)
    {
      sigma[,,i]=Sigma1[,,i]%*%(solve(V)[((i-1)*2*p+1):(i*2*p),((k-1)*2*p+1):(k*2*p)])
    }
    Vs[,,k]=apply(sigma,c(1,2),sum)
  }
  Vss=matrix(0,2*p,2*p)
  for (i in 1:K)
  {
    for (j in 1:K)
    {
      Vss=Vss+Sigma1[,,i]%*%(solve(V)[((i-1)*2*p+1):(i*2*p),((j-1)*2*p+1):(j*2*p)])%*%Sigma1[,,j]
    }
  }
  ##########wopt[,,k]=pi_{k,opt}#########
  wopt=array(0,dim=c(2*p,2*p,K))
  woptstar=array(0,dim=c(2*p,2*p,K))
  for (k in 1:K)
  {
    wopt[,,k]=solve(Vss)%*%Vs[,,k]%*%Sigma1[,,k]
    woptstar[,,k]=solve(Vss)%*%Vs[,,k]
  }
  ##########asd#########
  nvopt=solve(Vss)
  Beta_asd=sqrt((diag(nvopt/N))[1:p])  #asd of beta in vector
  Phi_asd=sqrt((diag(nvopt/N))[(p+1):(2*p)])  #asd of phi in vector

  return(list(wopt=wopt,woptstar=woptstar,Omega0=Omega0,Sigma1=Sigma1,nvopt=nvopt,Beta_asd=Beta_asd,Phi_asd=Phi_asd))
}
psi<-function(x,tau)
{
  tau-1*(x<0)
}
rk<-function(k,ehat,gehat,tau,N)
{
  #tau is a scalar; k is a lag
  btau=quantile(ehat,tau)[[1]]
  mug=mean(gehat)
  varg=var(gehat)
  rk=sum(psi(ehat[(k+1):length(ehat)]-btau,tau)*(gehat[(1):(length(gehat)-k)]-mug))/N/sqrt((tau-tau^2)*varg)
  return(rk)
}
OmegaKK<-function(X,XXt,wt,woptstar,Omega0,Sigma1,nvopt,beta,ehat,gehat,Tau,p,N,KK)
{
  #p is the dimension of beta; woptstar: array with dim=c(2*p,2*p,K)
  #nvopt: 2p*2p covariance matrix
  btau=as.vector(quantile(ehat,Tau))  #K*1 vector
  mug=mean(gehat)
  varg=apply(gehat,2,var)
  h=0.9*min(as.numeric(quantile(ehat,0.75)-quantile(ehat,0.25))/1.34,sd(ehat))*N^(-1/5)  #bandwidth
  density=density(btau,ehat,N,h,K)  #K*1 vector
  gamma=gamma(K,Tau,density)  #K*K matrix
  ##########GM#########
  GM=matrix(embed(gehat[-length(gehat)],KK),N-lagm-KK,KK) #(N-lagm-KK)*KK
  ##########DM#########
  d=array(0,dim=c(2*p,KK,N-lagm-KK,K))
  for (k in 1:K)
  {
    for (i in 1:(N-lagm-KK))  #t: from p+KK+1 to N
    {
      Xa=btau[k]/as.numeric(1+abs(X[i+KK,])%*%beta)*matrix(abs(X[i+KK,]),nrow=p,ncol=1)
      Xb=1/as.numeric(1+abs(X[i+KK,])%*%beta)*matrix(X[i+KK,],nrow=p,ncol=1)
      d[,,i,k]=matrix(rbind(Xa,Xb),nrow=2*p,ncol=1)%*%matrix(GM[i,]-mug,nrow=1,ncol=KK)
    }
  }
  DM=apply(d,c(1,2,4),mean)  #array with dim=c(2*p,KK,K)
  ##########Omega3#########
  D=matrix(0,2*p+1,KK)
  for (i in 1:(N-lagm-KK))
  {
    D=D+wt[i+KK]*t(matrix(XXt[i+KK,],nrow=1))%*%matrix(GM[i,]-mug,nrow=1,ncol=KK)
  }
  Omega3=D/N
  ###########Sigma3#########
  Sigma3=array(0,dim=c(2*p,2*p+1,K))
  sigma=array(0,dim=c(2*p,2*p+1,K))
  for (i in 1:K)
  {
    for (k in 1:K)
    {
      sigma[,,k]=gamma[i,k]*woptstar[,,k]%*%rbind(cbind(-beta,diag(p),matrix(0,p,p)),cbind(matrix(0,p,p+1),diag(p)))
    }
    Sigma3[,,i]=apply(sigma,c(1,2),sum)
  }
  ##########OmegaKK#########
  Omega=matrix(0,nrow=KK*K,ncol=KK*K)
  for (i in 1:K)
  {
    for (j in 1:K)
    {
      c=density[i]*density[j]/sqrt((Tau[i]-Tau[i]^2)*(Tau[j]-Tau[j]^2))/varg
      A1=t(DM[,,i])%*%Sigma3[,,j]%*%solve(Omega0)%*%Omega3
      A2=t(Omega3)%*%solve(Omega0)%*%t(Sigma3[,,i])%*%DM[,,j]
      B=t(DM[,,i])%*%nvopt%*%DM[,,j]
      Omega[((i-1)*KK+1):(i*KK),((j-1)*KK+1):(j*KK)]=c*(gamma[i,j]*varg*diag(KK)-A1-A2+B)
    }
  }
  return(Omega)
}

###########################################################################
# The following package is used to do order selection for LDAR(p) by BICs #
###########################################################################

ldarBIC<-function(N,Y,K,pmax)
{
  #####BIC based on self-weighted and doubly-weighted estimators####
  # N: sample size: Y: data to fit LDAR model                      #
  # K: the no. of quantiile levels                                 #
  # pmax: the maximum order in selection                           #
  ##################################################################
  BIC1<-BIC2<-BIC3<-NULL
  BIC4<-matrix(0,pmax,K)
  Tau=(1:K)/(K+1)
  #Create Responose and Regressor
  target<-Y[(pmax+1):N]         #Responose
  X=as.matrix((embed(Y[-length(Y)],pmax)))  #(N-pmax)*pmax
  XX=cbind(abs(X),X)      #Regressor without intercept (N-pmax)*2pmax
  XXt=cbind(1,XX)   #Regressor with intercept (N-pmax)*(2pmax+1)
  ###############LDAR(p) & self-weighted estimation#################
  w0=1/(1+apply(abs(X),1,sum)) #initial self-weights
  fit0=rq(target ~ XX, tau=Tau, weights=w0)
  est0=coef(fit0) #(2p+1)*K
  b_tilde1=matrix(est0[1,],1,K)
  beta_tilde1=matrix(est0[2:(pmax+1),],pmax,K)
  Phi_tilde1=matrix(est0[(pmax+2):(2*pmax+1),],pmax,K)
  Beta_tilde=apply(matrix(abs(beta_tilde1),pmax,K),1,sum)/sum(abs(b_tilde1))

  w1=as.vector(1/(1+abs(X)%*%(Beta_tilde+10^(-5))))
  w2=as.vector(1/(1+abs(X)%*%Beta_tilde)) #optimal self-weights
  fit1=rq(target ~ XX, tau=Tau, weights=w1)
  est=coef(fit1) #(2p+1)*K
  fit2=rq(target ~ XX, tau=Tau, weights=w2)

  b_tilde=matrix(est[1,],1,K)
  beta_tilde=matrix(est[2:(pmax+1),],pmax,K)
  Phi_tilde=matrix(est[(pmax+2):(2*pmax+1),],pmax,K)
  ########################doublely weighted estimation#####################
  #initial values for optimal weight
  Betahat0=apply(matrix(abs(beta_tilde),pmax,K),1,sum)/sum(abs(b_tilde))
  Phihat0=apply(matrix(Phi_tilde,pmax,K),1,mean)
  #optimal star weight
  opt<-asd(pmax,N,K,target,X,XXt,Betahat0,Phihat0,b_tilde,Tau,w1)
  woptstar=opt$woptstar
  #####################doubly weighted estimator based on optimal weight#################
  BetaPhi<-matrix(0,2*pmax,1)
  for (k in 1:K)
  {
    BetaPhi=BetaPhi+woptstar[,,k]%*%rbind(matrix(beta_tilde[,k],nrow=pmax),matrix(Phi_tilde[,k],nrow=pmax))
  }
  Betahat=BetaPhi[1:pmax,]  #pmax*1
  Phihat=BetaPhi[(pmax+1):(2*pmax),]  #pmax*1

  residual1=residuals(fit1)
  residual2=residuals(fit2)
  ehat=(target-X%*%Phihat)/(1+abs(X)%*%Betahat)  #(N-pmax)*1
  bhat=as.vector(quantile(ehat,Tau))  #K*1
  lambdaK=matrix(rep(BetaPhi,K),2*pmax,K)*rbind(t(matrix(rep(bhat,pmax),K,pmax)),matrix(1,pmax,K))
  residual3=matrix(rep(target,K),N-pmax,K)-XXt%*%rbind(bhat,lambdaK) #(N-pmax)*K; save e_t-Q_tau(e_t|F_{t-1})

  TauK=t(matrix(rep(Tau,N-pmax),K,N-pmax)) #(N-pmax)*K
  wthatK1=matrix(rep(w1,K),N-pmax,K) #(N-pmax)*K
  wthatK2=matrix(rep(w2,K),N-pmax,K) #(N-pmax)*K
  rhoresid1=residual1*(TauK-1*(residual1<0))*wthatK1 #(N-pmax)*K
  rhoresid2=residual2*(TauK-1*(residual2<0))*wthatK2 #(N-pmax)*K
  rhoresid3=residual3*(TauK-1*(residual3<0))*wthatK1 #(N-pmax)*K

  sigmahat1=apply(rhoresid1,2,mean) #K*1
  sigmahat2=apply(rhoresid2,2,mean) #K*1
  sigmahat3=apply(rhoresid3,2,mean) #K*1

  BICfit1=2*(N-pmax)*mean(log(sigmahat1))
  BICfit2=2*(N-pmax)*mean(log(sigmahat2))
  BICfit3=2*(N-pmax)*mean(log(sigmahat3))
  BICfit4=2*(N-pmax)*log(sigmahat1)

  pc=(2*pmax+1)*log(N-pmax)
  BIC1[pmax]=BICfit1+pc
  BIC2[pmax]=BICfit2+pc
  BIC3[pmax]=BICfit3+pc
  BIC4[pmax,]=BICfit4+pc

  for (p in 1:(pmax-1))
  {
    #Create Regressor
    X=matrix(embed(Y[-length(Y)],p)[-(1:(pmax-p)),],N-pmax,p)  #(N-pmax)*p
    XX=cbind(abs(X),X)    #Regressor without intercept (N-pmax)*2p
    XXt=cbind(1,XX)   #Regressor with intercept (N-pmax)*(2p+1)
    ######################LDAR(p) & self-weighted estimation#########################
    fit1=rq(target ~ XX, tau=Tau, weights=w1)
    est=coef(fit1) #(2p+1)*K
    fit2=rq(target ~ XX, tau=Tau, weights=w2)

    b_tilde=matrix(est[1,],1,K)
    beta_tilde=matrix(est[2:(p+1),],p,K)
    Phi_tilde=matrix(est[(p+2):(2*p+1),],p,K)
    ########################doublely weighted estimation#####################
    #initial values for optimal weight
    Betahat0=apply(matrix(abs(beta_tilde),p,K),1,sum)/sum(abs(b_tilde))
    Phihat0=apply(matrix(Phi_tilde,p,K),1,mean)
    #optimal star weight
    opt<-asd(p,N,K,target,X,XXt,Betahat0,Phihat0,b_tilde,Tau,w1)
    woptstar=opt$woptstar
    #####################doubly weighted estimator based on optimal weight#################
    BetaPhi<-matrix(0,2*p,1)
    for (k in 1:K)
    {
      BetaPhi=BetaPhi+woptstar[,,k]%*%rbind(matrix(beta_tilde[,k],nrow=p),matrix(Phi_tilde[,k],nrow=p))
    }
    Betahat=BetaPhi[1:p,]  #p*1
    Phihat=BetaPhi[(p+1):(2*p),]  #p*1

    residual1=residuals(fit1)
    residual2=residuals(fit2)
    ehat=(target-X%*%Phihat)/(1+abs(X)%*%Betahat)  #(N-pmax)*1
    bhat=as.vector(quantile(ehat,Tau))  #K*1
    lambdaK=matrix(rep(BetaPhi,K),2*p,K)*rbind(t(matrix(rep(bhat,p),K,p)),matrix(1,p,K))
    residual3=matrix(rep(target,K),N-pmax,K)-XXt%*%rbind(bhat,lambdaK) #(N-pmax)*K; save e_t-Q_tau(e_t|F_{t-1})

    rhoresid1=residual1*(TauK-1*(residual1<0))*wthatK1 #(N-pmax)*K
    rhoresid2=residual2*(TauK-1*(residual2<0))*wthatK2 #(N-pmax)*K
    rhoresid3=residual3*(TauK-1*(residual3<0))*wthatK1 #(N-pmax)*K

    sigmahat1=apply(rhoresid1,2,mean) #K*1
    sigmahat2=apply(rhoresid2,2,mean) #K*1
    sigmahat3=apply(rhoresid3,2,mean) #K*1

    BICfit1=2*(N-pmax)*mean(log(sigmahat1))
    BICfit2=2*(N-pmax)*mean(log(sigmahat2))
    BICfit3=2*(N-pmax)*mean(log(sigmahat3))
    BICfit4=2*(N-pmax)*log(sigmahat1)

    pc=(2*p+1)*log(N-pmax)
    BIC1[p]=BICfit1+pc
    BIC2[p]=BICfit2+pc
    BIC3[p]=BICfit3+pc
    BIC4[p,]=BICfit4+pc
  }

  Ind1=which.min(BIC1)
  Ind2=which.min(BIC2)
  Ind3=which.min(BIC3)
  med<-NULL
  for (k in 1:K)
  {
    med[k]=which.min(BIC4[,k])
  }

  Ind4=min(med)
  Ind5=max(med)

  return(list(BIC1=BIC1,BIC2=BIC2,BIC3=BIC3,BIC4=BIC4,
              Ind1=Ind1,Ind2=Ind2,Ind3=Ind3,Ind4=Ind4,Ind5=Ind5))
}

############################################################################
# The following package is used to do estimation and diagnosis for LDAR(p) #
############################################################################

ldarest=function(n,p,lagm,K,KK,ell,Y)
{
  ###############Estimation and do diagnostic checking for LDAR(p) model###############
  # n: sample size; p: order in LDAR model;                                           #
  # lagm: largest lag in model(for non-sparsity model,lagm=p);                        #
  # K: number of quantile levels in doubly weighted CQE;                              #
  # KK: maximum lag in QACF and Q_{BP} test;                                          #
  # ell: number of the simulation replication to obtain critical values;              #
  # Y: data to fit LDAR(p) model                                                      #
  #####################################################################################

  ##############################Data generation#############################
  #Create Responose and Regressor
  target=Y[(p+1):length(Y)]         #Responose
  X=as.matrix((embed(Y[-length(Y)],p)))
  XX=cbind(abs(X),X)     #Regressor without intercept
  XXt=cbind(rep(1,n-p),XX)

  #######################################Estimation################################
  Tau=(1:K)/(K+1)
  ###################LDAR(p) & self-weighted estimation#####################
  w1=1/(1+apply(abs(X),1,sum))
  fit1=rq(target ~ XX, tau=Tau, weights=w1)
  est1=coef(fit1)
  b_tilde1=est1[1,]
  beta_tilde1=est1[2:(p+1),]
  Phi_tilde1=est1[(p+2):(2*p+1),]

  Beta_tilde=apply(matrix(abs(beta_tilde1),p,length(Tau)),1,sum)/sum(abs(b_tilde1))
  w3=1/(1+abs(X)%*%abs(Beta_tilde))
  fit3=rq(target ~ XX, tau=Tau, weights=w3)
  est3=coef(fit3)
  b_tilde3=as.vector(est3[1,])
  beta_tilde3=matrix(est3[2:(p+1),],p,K)
  Phi_tilde3=matrix(est3[(p+2):(2*p+1),],p,K)
  ##########################double weighted estimation######################
  #initial values for optimal weight
  Betahat30=apply(matrix(abs(beta_tilde3),p,K),1,sum)/sum(abs(b_tilde3))
  Phihat30=apply(matrix(Phi_tilde3,p,K),1,mean)
  #optimal star weight
  w3_asd=asd(p,n,K,target,X,XXt,Betahat30,Phihat30,b_tilde3,Tau,w3)
  woptstar3=w3_asd$woptstar
  ###############doubly weighted estimation#################
  C3=matrix(0,2*p,1)
  for (k in 1:length(Tau))
  {
    C3=C3+woptstar3[,,k]%*%rbind(matrix(beta_tilde3[,k],nrow=p),matrix(Phi_tilde3[,k],nrow=p))
  }
  Phihat3=C3[(p+1):(2*p),]
  Betahat3=C3[1:p,]
  #########################obtain asd############################
  #optimal asd and covariance matrix nvopt
  w3_asdopt=asd(p,n,K,target,X,XXt,Betahat3,Phihat3,b_tilde3,Tau,w3)
  asd3Betaopt=w3_asdopt$Beta_asd
  asd3Phiopt=w3_asdopt$Phi_asd

  nvopt3=w3_asdopt$nvopt
  Omega03=w3_asdopt$Omega0
  Sigma13=w3_asdopt$Sigma1
  #########################diagnostic checking############################
  residual=(target-X%*%Phihat3)/(1+abs(X)%*%Betahat3)  #residual
  gCresid1=atan(residual)/pi+0.5
  gCresid2=atan(residual^2)/pi+0.5

  phoKC1=matrix(0,KK,length(Tau))
  phoKC2=matrix(0,KK,length(Tau))

  for (k in 1:K)
  {
    for (i in 1:KK)
    {
      phoKC1[i,k]=rk(i,residual,gCresid1,Tau[k],n)   #Cauchy
      phoKC2[i,k]=rk(i,residual,gCresid2,Tau[k],n)   #Cauchy
    }
  }
  RC1=apply(abs(phoKC1),1,max)
  RC2=apply(abs(phoKC2),1,max)

  QBPC1=n*sum(RC1^2)
  QBPC2=n*sum(RC2^2)

  #simulation for critical values or p-values of Q_BP
  SigmaC1=OmegaKK(X,XXt,w3,woptstar3,Omega03,Sigma13,nvopt3,Betahat3,residual,gCresid1,Tau,p,n,KK)
  SigmaC2=OmegaKK(X,XXt,w3,woptstar3,Omega03,Sigma13,nvopt3,Betahat3,residual,gCresid2,Tau,p,n,KK)

  SigmaC1PD=nearPD(SigmaC1)$mat
  SigmaC2PD=nearPD(SigmaC2)$mat

  #BTau: ell*(KK*K) matrix; each row is a sample
  Mu=rep(0,KK*length(Tau)) #mean
  set.seed(12345)
  BTauC1=mvrnorm(n = ell, Mu, SigmaC1PD, tol = 1e-6, empirical = FALSE)
  BTauC2=mvrnorm(n = ell, Mu, SigmaC2PD, tol = 1e-6, empirical = FALSE)
  BMC1=array(0,dim=c(KK,length(Tau),ell))
  BMC2=array(0,dim=c(KK,length(Tau),ell))
  for (i in 1:ell)
  {
    BMC1[,,i]=matrix(BTauC1[i,],KK,length(Tau))
    BMC2[,,i]=matrix(BTauC2[i,],KK,length(Tau))
  }
  asdRC1=apply(apply(abs(BMC1),c(1,3),max),1,sd)/sqrt(n)
  asdRC2=apply(apply(abs(BMC2),c(1,3),max),1,sd)/sqrt(n)
  #95% CI for suprho
  criticalRC1=apply(apply(abs(BMC1),c(1,3),max),1,quantile,probs=0.95)/sqrt(n)
  criticalRC2=apply(apply(abs(BMC2),c(1,3),max),1,quantile,probs=0.95)/sqrt(n)

  criticalC1=quantile(apply(apply(BMC1^2,c(1,3),max),2,sum),0.95)[[1]]
  criticalC2=quantile(apply(apply(BMC2^2,c(1,3),max),2,sum),0.95)[[1]]

  pvalue1=mean(apply(apply(BMC1^2,c(1,3),max),2,sum)>QBPC1)
  pvalue2=mean(apply(apply(BMC2^2,c(1,3),max),2,sum)>QBPC2)

  return(list(target=target,X=X,XX=XX,XXt=XXt,pvalue1=pvalue1,pvalue2=pvalue2,
              Betahat3=Betahat3,Phihat3=Phihat3,asd3Betaopt=asd3Betaopt,asd3Phiopt=asd3Phiopt,residual=residual,
              QBPC1=QBPC1,criticalC1=criticalC1,QBPC2=QBPC2,criticalC2=criticalC2,
              RC1=RC1,criticalRC1=criticalRC1,RC2=RC2,criticalRC2=criticalRC2,
              nvopt3=nvopt3))
}
