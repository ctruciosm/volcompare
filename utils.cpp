#include<RcppArmadillo.h>
#include<Rmath.h>

using namespace Rcpp;

// [[Rcpp::export]]
double dcsn_like(NumericVector params, NumericVector r){
  int n = r.size();
  NumericVector lambda(n), lambda_star(n), log_lik(n - 1);
  lambda_star[0] = 0.0;
  lambda[0] = params[0];
  for(int i = 1; i < n; i++){
    lambda_star[i] = params[1] * lambda_star[i - 1] + params[2] * (pow(r[i - 1], 2) / exp(2 * lambda[i - 1]) - 1);
    lambda[i] = params[0] + lambda_star[i];
    log_lik[i - 1] = 0.5* pow(r[i], 2) / exp(2 * lambda[i]) + lambda[i];
  }
  return sum(log_lik);
}

inline int sign(double x){
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
double dcsn_lev_like(NumericVector params, NumericVector r){
  int n = r.size();
  double u;
  NumericVector lambda(n), lambda_star(n), log_lik(n - 1);
  lambda_star[0] = 0.0;
  lambda[0] = params[0];
  double kapa = exp(params[2]);
  double kapastar = kapa / (1.0 + exp(-params[3]));
  for(int i = 1; i < n; i++){
    u = pow(r[i - 1], 2) / exp(2 * lambda[i - 1]) - 1;
    lambda_star[i] = params[1] * lambda_star[i - 1] + kapa * u + kapastar * sign(-r[i - 1]) * (u + 1);
    lambda[i] = params[0] + lambda_star[i];
    log_lik[i - 1] = 0.5* pow(r[i], 2) / exp(2 * lambda[i]) + lambda[i];
  }
  return sum(log_lik);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_dcsn(NumericVector y){
  NumericVector coeff(3), vi(3);
  double omega, phi, kapa;
  double omega_min = -0.01, omega_max= 0.05, phi_min = 0.25,  phi_max = 0.99, kapa_min = 0.05,  kapa_max = 0.25, n_omega = 5, n_phi = 5,  n_kapa = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min) / n_omega;
  double lm_phi = (phi_max - phi_min) / n_phi;
  double lm_kapa = (kapa_max - kapa_min) / n_kapa;

  for(int no = 0; no < n_omega; no++){
    for(int np = 0; np < n_phi; np++){
      for(int nk =0; nk < n_kapa; nk++){
        omega = omega_min + no * lm_omega;
        phi = phi_min + np * lm_phi; 
        kapa = kapa_min + nk * lm_kapa;
        coeff[0] = omega;
        coeff[1] = phi;
        coeff[2] = kapa;
        nml = dcsn_like(coeff, y);
        if (nml < ml){
          vi[0] = coeff[0];
          vi[1] = coeff[1];
          vi[2] = coeff[2];
          ml=nml;
        }
      }
    }
  }
  return(vi);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_dcsn_lev(NumericVector y){
  NumericVector coeff(4), vi(4);
  double omega, phi, kapa, kapa2;
  double omega_min = -0.01, omega_max= 0.05, phi_min = 0.25,  phi_max = 0.99, kapa_min = 0.05,  kapa_max = 0.25, kapa2_min = 0.01, kapa2_max = 0.20, n_omega = 5, n_phi = 5,  n_kapa = 5, n_kapa2 = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min) / n_omega;
  double lm_phi = (phi_max - phi_min) / n_phi;
  double lm_kapa = (kapa_max - kapa_min) / n_kapa;
  double lm_kapa2 = (kapa2_max - kapa2_min) / n_kapa2;
  
  for(int no = 0; no < n_omega; no++){
    for(int np = 0; np < n_phi; np++){
      for(int nk =0; nk < n_kapa; nk++){
        for(int nk2 =0; nk2 < n_kapa2; nk2++){
          omega = omega_min + no * lm_omega;
          phi = phi_min + np * lm_phi; 
          kapa = kapa_min + nk * lm_kapa;
          kapa2 = kapa2_min + nk2 * lm_kapa2;
          coeff[0] = omega;
          coeff[1] = phi;
          coeff[2] = log(kapa);
          coeff[3] = log(kapa2 / (kapa - kapa2));
          nml = dcsn_lev_like(coeff, y);
          if (nml < ml){
            vi[0] = coeff[0];
            vi[1] = coeff[1];
            vi[2] = coeff[2]; 
            vi[3] = coeff[3]; 
            ml=nml;
          }
        }
      }
    }
  }
  return(vi);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]  
arma::mat QL(arma::mat VaR, arma::vec r, double alpha){
  int n = VaR.n_rows, k = VaR.n_cols;
  arma::mat value(n,k);
  for (int j=0; j<n; j++){
    for (int i=0; i<k; i++){
      if(r[j]<=VaR(j,i)){
        value(j,i) = (alpha-1)*(r[j]-VaR(j,i));
      } else{
        value(j,i) = alpha*(r[j]-VaR(j,i));
      }
    }
  }
  return(value);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]  
arma::mat AL(arma::mat VaR, arma::mat ES, arma::vec r, double alpha){
  int n = VaR.n_rows, k = VaR.n_cols;
  arma::mat value(n,k);
  for (int j=0; j<n; j++){
    for (int i=0; i<k; i++){
      if(r[j]<=VaR(j,i)){
        value(j,i) =((-1/ES(j,i))*(ES(j,i) - VaR(j,i) + (VaR(j,i) - r[j])/alpha) - (-log(-ES(j,i))) + (1-log(1-alpha)));
      } else{
        value(j,i) =((-1/ES(j,i))*(ES(j,i) - VaR(j,i)) - (-log(-ES(j,i))) + (1-log(1-alpha)));
      }
    }
  }
  return(value);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]  
arma::mat NZ(arma::mat VaR, arma::mat ES, arma::vec r, double alpha){
  int n = VaR.n_rows, k = VaR.n_cols;
  arma::mat value(n,k);
  for (int j=0; j<n; j++){
    for (int i=0; i<k; i++){
      if(r[j]<=VaR(j,i)){
        value(j,i) = ((1/(2*sqrt(-ES(j,i))))*(ES(j,i) - VaR(j,i) + (VaR(j,i) - r[j])/alpha) + sqrt(-ES(j,i)));
      } else{
        value(j,i) = ((1/(2*sqrt(-ES(j,i))))*(ES(j,i) - VaR(j,i)) + sqrt(-ES(j,i)));
      }
    }
  }
  return(value);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]  
arma::mat FZG(arma::mat VaR, arma::mat ES, arma::vec r, double alpha){
  int n = VaR.n_rows, k = VaR.n_cols;
  arma::mat value(n,k);
  for (int j=0; j<n; j++){
    for (int i=0; i<k; i++){
      if(r[j]<=VaR(j,i)){
        value(j,i) = (1 - alpha)*VaR(j,i) - r[j] + (exp(ES(j,i))/(1+exp(ES(j,i))))*(ES(j,i) - VaR(j,i) + (VaR(j,i) - r[j])/alpha) - log(1+exp(ES(j,i))) + log(2);
      } else{
        value(j,i) = -alpha*VaR(j,i) + (exp(ES(j,i))/(1+exp(ES(j,i))))*(ES(j,i) - VaR(j,i)) - log(1+exp(ES(j,i))) + log(2);
      }
    }
  }
  return(value);
}
