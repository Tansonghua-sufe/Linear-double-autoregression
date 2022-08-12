// [[Rcpp::depends(RcppEigen)]]
#include<RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::export]]
Eigen::VectorXd segment_e(Eigen::VectorXd x,int location,int num){
  // ======================
  // auxiliary function: return x[location,location+1,...,location+num-1]
  // ======================
  Eigen::VectorXd y;
  Eigen::VectorXd z;
  z.resize(num);
  y=x.segment(location,num);
  for(int i=0;i<num;i=i+1){
    z(i)=y(num-i-1);
  }
  return z;
}
Eigen::VectorXd rep(double value,int number){
  // ======================
  // auxiliary function: replicate elements
  // ======================
  Eigen::VectorXd y;
  y.resize(number);y.fill(value);
  return y;
}



//[[Rcpp::export]]
Eigen::VectorXd DGP (Eigen::VectorXd epsilon,Eigen::VectorXd lambda_1,Eigen::VectorXd lambda_2) {
  // ===============================================
  // Data Generation Process
  // epsilon: innovations in model
  // lambda_1: conditional location parameter
  // lambda_2: conditional volatility parameter
  // ===============================================
  double p = lambda_1.size();
  double n = epsilon.size();
  double mean;
  double varepsilon;
  Eigen::VectorXd Y = epsilon;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;

  Y_2.resize(p+1); Y_2.fill(0);
  (Y.head(p)).fill(0);

    for(int t=p; t<n;t=t+1)
    {
      Y_1 = segment_e(Y,t-p,p);
      Y_2(0)= 1;
      Y_2.segment(1,p)= Y_1.cwiseAbs();
      mean=lambda_1.adjoint()*Y_1;
      varepsilon=lambda_2.adjoint()*Y_2;
      Y(t) = mean + epsilon(t)*varepsilon;
    }
  return Y;
}

//[[Rcpp::export]]
double likehood_LDAR_EQMLE(Eigen::VectorXd lambda,Eigen::VectorXd Y) {
  // ==========================================
  // Negative log-likelihood function of EMQLE
  // lambda is coefficient
  // Y is observation
  // ==========================================
  double p = (lambda.size()-1)/2;
  double n = Y.size();
  double likehood = 0;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1= lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);
  Y_2.resize(p+1); Y_2.fill(0);
  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();
    likehood = likehood+log(lambda_2.adjoint()*Y_2)+fabs((Y(t)-lambda_1.adjoint()*Y_1))/
        (lambda_2.adjoint()*Y_2);
  }
  return likehood;
}

//[[Rcpp::export]]
double likehood_LDAR_GQMLE(Eigen::VectorXd lambda,Eigen::VectorXd Y) {
  // ==========================================
  // Negative log-likelihood function of GMQLE
  // lambda is coefficient
  // Y is observation
  // ==========================================
  double p = (lambda.size()-1)/2;
  double n = Y.size();
  double likehood = 0;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1= lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);
  Y_2.resize(p+1); Y_2.fill(0);
  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();
    likehood = likehood+log(lambda_2.adjoint()*Y_2)+pow((Y(t)-lambda_1.adjoint()*Y_1),2)/
        (2*pow((lambda_2.adjoint()*Y_2),2));
  }
  return likehood;
}


//[[Rcpp::export]]
Eigen::VectorXd resdual(Eigen::VectorXd Y,Eigen::VectorXd fitted_lambda) {
  // ========================================================
  // Extract residuals of LDAR model based on the given parameter
  // Y is observation
  // fitted_lambda is the given parameter
  // ========================================================
  double p = (fitted_lambda.size()-1)/2;
  double n = Y.size();

  Eigen::VectorXd residuals=Y;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1;
  Eigen::VectorXd lambda_2;

  lambda_1=fitted_lambda.segment(0,p);
  lambda_2=fitted_lambda.segment(p,p+1);
  Y_2.resize(p+1); Y_2.fill(0);
  (residuals.head(p)).fill(0);

  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0) = 1;Y_2.segment(1,p)= Y_1.cwiseAbs();
    residuals(t) = (Y(t)-lambda_1.adjoint()*Y_1)/(lambda_2.adjoint()*Y_2);
  }
  return residuals;
}

//[[Rcpp::export]]
Eigen::MatrixXd Omega_EQMLE(Eigen::VectorXd lambda_1,Eigen::VectorXd lambda_2,Eigen::VectorXd Y,double E_1,double E_2) {
  // =========================================================
  // Extract Omega matrix of covariance matrix based on EQMLE
  // lambda_1 and lambda_2 are the estimation of EQMLE
  // E_1 is the expectation of the residuals
  // E_2 is the expected variance of the residuals
  // =========================================================
  const double p = lambda_1.size();
  const double n = Y.size();
  Eigen::MatrixXd omega;
  Eigen::MatrixXd omega_11;
  Eigen::MatrixXd omega_12;
  Eigen::MatrixXd omega_22;

  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;

  Y_2.resize(p+1); Y_2.fill(0);
  omega.resize(2*p+1,2*p+1);
  omega_11.resize(p,p);omega_12.resize(p,p+1);omega_22.resize(p+1,p+1);
  omega_11.fill(0);omega_12.fill(0);omega_22.fill(0);

  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();
    omega_11 = omega_11+Y_1*Y_1.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
    omega_12 = omega_12+E_1*Y_1*Y_2.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
    omega_22 = omega_22+(E_2-1)*Y_2*Y_2.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
  }

  omega.topLeftCorner(p,p)=omega_11;omega.topRightCorner(p,p+1)=omega_12;
  omega.bottomLeftCorner(p+1,p)=omega_12.adjoint();omega.bottomRightCorner(p+1,p+1)=omega_22;

  omega = omega.array()/n;

  return omega;
}

//[[Rcpp::export]]
Eigen::MatrixXd Sigma_EQMLE(Eigen::VectorXd lambda,Eigen::VectorXd Y,double f_0) {
  // ============================================================
  // Extract Sigma matrix of covariance matrix based on EQMLE
  // lambda is the estimation of EQMLE
  // f_0 is the density function of innovation at zero, i.e. f(0)
  // ============================================================
  const double p = (lambda.size()-1)/2;
  const double n = Y.size();
  Eigen::MatrixXd sigma;
  Eigen::MatrixXd s_11;
  Eigen::MatrixXd s_12;
  Eigen::MatrixXd s_22;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1=lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);

  sigma.resize(2*p+1,2*p+1);
  s_11.resize(p,p);s_12.resize(p,p+1);s_22.resize(p+1,p+1);

  s_11.fill(0);s_12.fill(0);s_22.fill(0);

  Y_2.resize(p+1); Y_2.fill(0);
  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();

    s_11 = s_11+f_0*Y_1*Y_1.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
    s_22 = s_22+0.5*Y_2*Y_2.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
  }

  sigma.topLeftCorner(p,p)=s_11;sigma.topRightCorner(p,p+1)=s_12;
  sigma.bottomLeftCorner(p+1,p)=s_12.adjoint();sigma.bottomRightCorner(p+1,p+1)=s_22;

  sigma = sigma.array()/n;

  return sigma;
}

//[[Rcpp::export]]
Eigen::MatrixXd Xi_EQMLE(Eigen::VectorXd lambda,Eigen::VectorXd Y,double E_1,double E_2,double f_0) {
  // ============================================================
  // Extract covariance matrix based on EQMLE
  // lambda is the estimation of EQMLE
  // E_1 is the expectation of the residuals
  // E_2 is the expected variance of the residuals
  // f_0 is the density function of innovation at zero, i.e. f(0)
  // ============================================================
  const double p = (lambda.size()-1)/2;
  Eigen::MatrixXd Xi;
  Eigen::MatrixXd sigma;
  Eigen::MatrixXd omega;
  Eigen::VectorXd lambda_1=lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);

  omega = Omega_EQMLE(lambda_1,lambda_2,Y,E_1,E_2);
  sigma = Sigma_EQMLE(lambda,Y,f_0);


  Xi.resize(2*p+1,2*p+1);
  Xi = (sigma.inverse())*omega*(sigma.inverse())/4;
  return Xi;
}


//[[Rcpp::export]]
Eigen::MatrixXd Omega_GQMLE(Eigen::VectorXd lambda_1,Eigen::VectorXd lambda_2,Eigen::VectorXd Y,double E_3,double E_4) {
  // =========================================================
  // Extract Omega matrix of covariance matrix based on GQMLE
  // lambda_1 and lambda_2 are the estimation of GQMLE
  // E_3 is the third moment of the residuals
  // E_4 is the fourth moment of the residuals - 1
  // =========================================================
  const double p = lambda_1.size();
  const double n = Y.size();
  Eigen::MatrixXd omega;
  Eigen::MatrixXd omega_11;
  Eigen::MatrixXd omega_12;
  Eigen::MatrixXd omega_22;

  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;

  Y_2.resize(p+1); Y_2.fill(0);
  omega.resize(2*p+1,2*p+1);
  omega_11.resize(p,p);omega_12.resize(p,p+1);omega_22.resize(p+1,p+1);
  omega_11.fill(0);omega_12.fill(0);omega_22.fill(0);

  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();
    omega_11 = omega_11+Y_1*Y_1.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
    omega_12 = omega_12+E_3*Y_1*Y_2.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
    omega_22 = omega_22+(E_4-1)*Y_2*Y_2.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
  }


  omega.topLeftCorner(p,p)=omega_11;omega.topRightCorner(p,p+1)=omega_12;
  omega.bottomLeftCorner(p+1,p)=omega_12.adjoint();omega.bottomRightCorner(p+1,p+1)=omega_22;

  omega = omega.array()/n;

  return omega;
}

//[[Rcpp::export]]
Eigen::MatrixXd Sigma_GQMLE(Eigen::VectorXd lambda,Eigen::VectorXd Y) {
  // ============================================================
  // Extract Sigma matrix of covariance matrix based on GQMLE
  // lambda is the estimation of GQMLE
  // ============================================================
  const double p = (lambda.size()-1)/2;
  const double n = Y.size();
  Eigen::MatrixXd hesse;
  Eigen::MatrixXd h_11;
  Eigen::MatrixXd h_12;
  Eigen::MatrixXd h_22;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1=lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);

  hesse.resize(2*p+1,2*p+1);
  h_11.resize(p,p);h_12.resize(p,p+1);h_22.resize(p+1,p+1);

  h_11.fill(0);h_12.fill(0);h_22.fill(0);


  Y_2.resize(p+1); Y_2.fill(0);
  for(int t=p; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();

    h_11 = h_11+Y_1*Y_1.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
    h_22 = h_22+2*Y_2*Y_2.adjoint()/pow((lambda_2.adjoint()*Y_2),2);
  }

  hesse.topLeftCorner(p,p)=h_11;hesse.topRightCorner(p,p+1)=h_12;
  hesse.bottomLeftCorner(p+1,p)=h_12.adjoint();hesse.bottomRightCorner(p+1,p+1)=h_22;

  hesse = hesse.array()/n;

  return hesse;
}

//[[Rcpp::export]]
Eigen::MatrixXd Xi_GQMLE(Eigen::VectorXd lambda,Eigen::VectorXd Y,double E_3,double E_4) {
  // ============================================================
  // Extract covariance matrix based on GQMLE
  // lambda is the estimation of GQMLE
  // E_3 is the third moment of the residuals
  // E_4 is the fourth moment of the residuals - 1
  // ============================================================
  const double p = (lambda.size()-1)/2;
  Eigen::MatrixXd Xi;
  Eigen::MatrixXd sigma;
  Eigen::MatrixXd omega;
  Eigen::VectorXd lambda_1=lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);
  sigma = Sigma_GQMLE(lambda,Y);
  omega = Omega_GQMLE(lambda_1,lambda_2,Y,E_3,E_4);

  Xi.resize(2*p+1,2*p+1);
  Xi = (sigma.inverse())*omega*(sigma.inverse());
  return Xi;
}


//[[Rcpp::export]]
double likehood_EQMLE_p_max(Eigen::VectorXd lambda,Eigen::VectorXd Y, int p_max) {
  // ========================================================
  // Negative log-likelihood function of EMQLE
  // (start from p_max and only used in BIC)
  // lambda is coefficient
  // Y is observation
  // p_max is the max lag, predetermined positive integer.
  // ========================================================
  double p = (lambda.size()-1)/2;
  double n = Y.size();
  double likehood = 0;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1= lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);
  Y_2.resize(p+1); Y_2.fill(0);

  for(int t=p_max; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();
    likehood = likehood+log(lambda_2.adjoint()*Y_2)+fabs((Y(t)-lambda_1.adjoint()*Y_1))/
        (lambda_2.adjoint()*Y_2);
  }

  return likehood;
}

//[[Rcpp::export]]
double likehood_LDAR_GQMLE_p_max(Eigen::VectorXd lambda,Eigen::VectorXd Y, int p_max) {
  // ========================================================
  // Negative log-likelihood function of GMQLE
  // (start from p_max and only used in BIC)
  // lambda is coefficient
  // Y is observation
  // p_max is the max lag, predetermined positive integer.
  // ========================================================
  double p = (lambda.size()-1)/2;
  double n = Y.size();
  double likehood = 0;
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1= lambda.segment(0,p);
  Eigen::VectorXd lambda_2=lambda.segment(p,p+1);
  Y_2.resize(p+1); Y_2.fill(0);

  for(int t=p_max; t<n;t=t+1)
  {
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0)= 1;
    Y_2.segment(1,p)= Y_1.cwiseAbs();
    likehood = likehood+log(lambda_2.adjoint()*Y_2)+pow((Y(t)-lambda_1.adjoint()*Y_1),2)/
        (2*pow((lambda_2.adjoint()*Y_2),2));
  }

  return likehood;
}

//[[Rcpp::export]]
double indicate(double eta){
  // ==============================================
  // auxiliary function used in portmanteau test
  // Indicative function
  // ==============================================
  double sgn_eta;
  if(eta>0){
    sgn_eta=1;
  }else{
    sgn_eta=-1;
  }
  return sgn_eta;
}

double cov(Eigen::VectorXd x,int k,int p){
  // ==============================================
  // auxiliary function used in portmanteau test
  // autocorrelation function
  //y = \sum_{t=p+k+1}(x_t-\bar{x})(x_{t-k}-\bar{x})
  // ==============================================
  double x_bar;
  int n=x.size();
  double y=0;
  x_bar=x.sum()/(n-p);
  for(int t=k+p;t<n;t=t+1){
    y=y+(x(t)-x_bar)*(x(t-k)-x_bar);
  }
  return y;
}

//[[Rcpp::export]]
Eigen::VectorXd sgn(Eigen::VectorXd x){
  // ==============================================
  // auxiliary function used in portmanteau test
  // sign function
  // ==============================================
  int num=x.size();
  Eigen::VectorXd y;
  y.resize(num);
  for(int i=0;i<num;i++){
    if(x(i)>0){
      y(i)=1.0;
    }else{
      if(x(i)==0){
        y(i)=0.0;
      }else{
        y(i)=-1.0;
      }
    }
  }
  return y;
}

Eigen::VectorXd pow_vector(Eigen::VectorXd x,int k){
  // ==============================================
  // auxiliary function used in portmanteau test
  // return x^k where x is vector
  // ==============================================
  Eigen::VectorXd y;
  int n=x.size();
  y.resize(n);
  for(int t=0;t<n;t=t+1){
    y(t)=pow(x(t),k);
  }
  return y;
}

//[[Rcpp::export]]
List portmanteau_test_EQMLE(Eigen::VectorXd fitted_lambda,Eigen::VectorXd Y,int m, Eigen::MatrixXd Sigma){
  // =========================================================
  // Portmanteau test of EQMLE
  // fitted_lambda is the estimation of EQMLE
  // Y is observation
  // m is predetermined positive integer
  // sigma is the Sigma matrix in covariance matrix of EQMLE
  // =========================================================
  int p = (fitted_lambda.size()-1)/2;
  int n = Y.size();
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1=fitted_lambda.segment(0,p);
  Eigen::VectorXd lambda_2=fitted_lambda.segment(p,p+1);
  Eigen::VectorXd res;
  Eigen::VectorXd res_abs;
  Eigen::VectorXd res_sgn;
  Eigen::VectorXd cf;
  Eigen::VectorXd v_t;
  Eigen::VectorXd G_t;
  double stat;
  double tau;
  double sigma_1;
  double sigma_2;
  double h;
  double eta_t;
  Eigen::MatrixXd V;
  Eigen::MatrixXd V_1;
  Eigen::MatrixXd V_2;
  Eigen::MatrixXd G;
  Eigen::MatrixXd Sigma_cf;

  Y_2.resize(p+1);
  cf.resize(2*m);cf.fill(0);
  v_t.resize(2*m+2*p+1);
  v_t.fill(0);
  G.resize(2*m+2*p+1,2*m+2*p+1);G.fill(0);

  G_t.resize(2*p+1);

  V_1.setIdentity(2*m,2*m);
  V_2.resize(2*m,2*p+1);V_2.fill(0);
  V.resize(2*m,2*p+1+2*m);

  Sigma_cf.resize(2*m,2*m);
  Sigma.resize(2*p+1,2*p+1);

  res=resdual(Y,fitted_lambda);
  res_sgn=sgn(res);
  res_abs=res.cwiseAbs();
  tau=res.sum()/(n-p);
  sigma_1=pow_vector(res,2).sum()/(n-p)-pow(tau,2);
  sigma_2=pow_vector(res,2).sum()/(n-p)-1;

  for(int t=p+m;t<n;t++){
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0) = 1;Y_2.segment(1,p) = Y_1.cwiseAbs();
    h = lambda_2.transpose()*Y_2;
    eta_t = res(t);

    G_t.segment(0,p)=-Y_1/h*indicate(eta_t);
    G_t.segment(p,p+1)=Y_2/h*(1-fabs(eta_t));
    v_t.segment(2*m,2*p+1)=-0.5*Sigma.inverse()*G_t;
    for(int k=0;k<m;k++){
      V_2.block(k,0,1,p)=V_2.block(k,0,1,p)-Y_1.transpose()*(res(t-k-1)-tau)/(h*sigma_1);
      V_2.block(k,p,1,p+1)=V_2.block(k,p,1,p+1)-tau*Y_2.transpose()*(res(t-k-1)-tau)/(h*sigma_1);
      V_2.block(k+m,p,1,p+1)=V_2.block(k+m,p,1,p+1)-(res_abs(t-k-1)-1)*Y_2.transpose()/(h*sigma_2);
      v_t(k)=(res(t)-tau)*(res(t-k-1)-tau)/sigma_1;
      v_t(k+m)=(fabs(res(t))-1)*(fabs(res(t-k-1))-1)/sigma_2;
    }

    G+=v_t*v_t.transpose();
  }

  V_2=V_2.array()/(n-p-m);
  V.leftCols(2*m)=V_1;V.rightCols(2*p+1)=V_2;

  G=G.array()/(n-p-m);

  Sigma_cf=V*G*V.transpose();

  for(int k=0;k<m;k++){
    cf(k)=cov(res,k+1,p)/cov(res,0,p);
    cf(k+m)=cov(res_abs,k+1,p)/cov(res_abs,0,p);
  }

  stat=n*cf.transpose()*(Sigma_cf).inverse()*cf;


  return(
    List::create(
      Named("rho")=cf.segment(0,m),         // autocorrelation coefficient of residuals
      Named("gamma")=cf.segment(m,m),       // autocorrelation coefficient of abs residuals
      Named("stat")=stat,                   // test statistic
      Named("Sigma_cf")=Sigma_cf           // covariance matrix, i.e. VGV^prime
    ));
}


//[[Rcpp::export]]
List portmanteau_test_GQMLE(Eigen::VectorXd fitted_lambda,Eigen::VectorXd Y,int m, Eigen::MatrixXd Sigma){
  // =========================================================
  // Portmanteau test of EQMLE
  // fitted_lambda is the estimation of EQMLE
  // Y is observation
  // m is predetermined positive integer
  // sigma is the Sigma matrix in covariance matrix of EQMLE
  // =========================================================
  int p = (fitted_lambda.size()-1)/2;
  int n = Y.size();
  Eigen::VectorXd Y_1;
  Eigen::VectorXd Y_2;
  Eigen::VectorXd lambda_1=fitted_lambda.segment(0,p);
  Eigen::VectorXd lambda_2=fitted_lambda.segment(p,p+1);
  Eigen::VectorXd res; // residuals
  Eigen::VectorXd res_abs; // absolute values of residuals
  Eigen::VectorXd res_sgn; // signs of resiudals
  Eigen::VectorXd cf;                 // autocorrelation coefficient of residuals and its absolute values
  Eigen::VectorXd v_t;                // need to G=E(v_t*v_t^\prime)
  Eigen::VectorXd G_t;
  double stat;                        // test statistic
  double tau_1;                       //\tau_1=E(sgn(\eta_t))
  double tau_2;                       //\tau_2=E(|\eta_t|)
  double sigma_xi;                    //\sigma_xi=1-(E(abs(eta_t)))^2
  double h;
  double eta_t;
  Eigen::MatrixXd V;
  Eigen::MatrixXd V_1;
  Eigen::MatrixXd V_2;
  Eigen::MatrixXd G;
  Eigen::MatrixXd Sigma_cf;     // covariance matrix, i.e. VGV^prime

  Y_2.resize(p+1);
  cf.resize(2*m);cf.fill(0);
  v_t.resize(2*m+2*p+1);
  v_t.fill(0);
  G.resize(2*m+2*p+1,2*m+2*p+1);G.fill(0);

  G_t.resize(2*p+1);

  V_1.setIdentity(2*m,2*m);
  V_2.resize(2*m,2*p+1);V_2.fill(0);
  V.resize(2*m,2*p+1+2*m);

  Sigma_cf.resize(2*m,2*m);
  Sigma.resize(2*p+1,2*p+1);

  res=resdual(Y,fitted_lambda);
  res_sgn=sgn(res);
  res_abs=res.cwiseAbs();
  tau_1=res_sgn.sum()/(n-p);
  tau_2=res_abs.sum()/(n-p);
  sigma_xi=1-pow(tau_2,2);
  Eigen::MatrixXd Sigma_inv = Sigma.inverse();Sigma_inv.resize(2*p+1,2*p+1);

  for(int t=p+m;t<n;t++){
    Y_1 = segment_e(Y,t-p,p);
    Y_2(0) = 1;Y_2.segment(1,p) = Y_1.cwiseAbs();
    h = lambda_2.transpose()*Y_2;
    eta_t = res(t);
    G_t.segment(0,p)=-Y_1/h*eta_t;
    G_t.segment(p,p+1)=Y_2/h*(1-pow(eta_t,2));
    v_t.segment(2*m,2*p+1)=-Sigma_inv*G_t;
    for(int k=0;k<m;k++){
      V_2.block(k,0,1,p)=V_2.block(k,0,1,p)-Y_1.transpose()*(res(t-k-1))/h;
      V_2.block(k+m,0,1,p)=V_2.block(k+m,0,1,p)-tau_1*Y_1.transpose()*(res_abs(t-k-1)-tau_2)/(h*sigma_xi);
      V_2.block(k+m,p,1,p+1)=V_2.block(k+m,p,1,p+1)-tau_2*(res_abs(t-k-1)-tau_2)*Y_2.transpose()/(h*sigma_xi);
      v_t(k)=res(t)*res(t-k-1); // the first part of v_t
      v_t(k+m)=(res_abs(t)-tau_2)*(res_abs(t-k-1)-tau_2)/sigma_xi;// the second part of v_t
    }

    G+=v_t*v_t.transpose();
  }

  V_2=V_2.array()/(n-p-m);
  V.leftCols(2*m)=V_1;V.rightCols(2*p+1)=V_2;

  G=G.array()/(n-p-m);

  Sigma_cf=V*G*V.transpose();

  for(int k=0;k<m;k++){
    cf(k)=cov(res,k+1,p)/cov(res,0,p);
    cf(k+m)=cov(res_abs,k+1,p)/cov(res_abs,0,p);
  }

  stat=n*cf.transpose()*(Sigma_cf).inverse()*cf;

  return(
    List::create(
      Named("rho")=cf.segment(0,m),         // autocorrelation coefficient of residuals
      Named("gamma")=cf.segment(m,m),       // autocorrelation coefficient of abs residuals
      Named("stat")=stat,                   // test statistic
      Named("Sigma_cf")=Sigma_cf           // covariance matrix, i.e. VGV^prime
    ));
}
