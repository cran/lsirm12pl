// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::export]]

Rcpp::List log_likelihood_cpp(arma::mat data, const int ndim, arma::mat beta_est, arma::mat theta_est, const double gamma_est, arma::mat z_est, arma::mat w_est,
                          const double missing){
  double log_likelihood = 0, dist_temp = 0;
  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  int i, j, k;
  arma::dmat dist(nsample,nitem,fill::zeros);

  dist.fill(0.0);
  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      dist_temp = 0.0;
      for(j = 0; j < ndim; j++) dist_temp += std::pow((z_est(k,j)-w_est(i,j)), 2.0);
      dist(k,i) = std::sqrt(dist_temp);
    }
  }


  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      if(data(k,i) != missing){
        if(data(k,i) == 1.0) log_likelihood += -std::log(1.0 + std::exp(-(beta_est(i) + theta_est(k) - gamma_est * dist(k,i))));
        else log_likelihood += -std::log(1.0 + std::exp(beta_est(i) + theta_est(k) - gamma_est * dist(k,i)));
      }
    }
  }


  Rcpp::List output;
  output["log_likelihood"] = log_likelihood;

  return(output);

} // function end

// [[Rcpp::export]]
Rcpp::List log_likelihood_normal_cpp(arma::mat data, const int ndim, arma::mat beta_est, arma::mat theta_est, const double gamma_est, arma::mat z_est, arma::mat w_est, const double sigma_est,
                              const double missing){
  double log_likelihood = 0, dist_temp = 0;
  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  int i, j, k;
  arma::dmat dist(nsample,nitem,fill::zeros);

  dist.fill(0.0);
  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      dist_temp = 0.0;
      for(j = 0; j < ndim; j++) dist_temp += std::pow((z_est(k,j)-w_est(i,j)), 2.0);
      dist(k,i) = std::sqrt(dist_temp);
    }
  }


  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      if(data(k,i) != missing){
        log_likelihood += -pow((data(k, i) - beta_est(i) - theta_est(k) + gamma_est * dist(k, i)), 2) / (2 * pow(sigma_est, 2));
      }
    }
  }


  Rcpp::List output;
  output["log_likelihood"] = log_likelihood;

  return(output);

} // function end



// [[Rcpp::export]]
Rcpp::List log_likelihood_normal2pl_cpp(arma::mat data, const int ndim, arma::mat beta_est, arma::mat alpha_est, arma::mat theta_est, const double gamma_est, arma::mat z_est, arma::mat w_est, const double sigma_est,
                                     const double missing){
  double log_likelihood = 0, dist_temp = 0;
  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  int i, j, k;
  arma::dmat dist(nsample,nitem,fill::zeros);

  dist.fill(0.0);
  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      dist_temp = 0.0;
      for(j = 0; j < ndim; j++) dist_temp += std::pow((z_est(k,j)-w_est(i,j)), 2.0);
      dist(k,i) = std::sqrt(dist_temp);
    }
  }


  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      if(data(k,i) != missing){
        log_likelihood += -pow((data(k, i) - beta_est(i) - alpha_est(i) * theta_est(k) + gamma_est * dist(k, i)), 2) / (2 * pow(sigma_est, 2));
      }
    }
  }


  Rcpp::List output;
  output["log_likelihood"] = log_likelihood;

  return(output);

} // function end

// [[Rcpp::export]]
Rcpp::List log_likelihood_2pl_cpp(arma::mat data, const int ndim, arma::mat beta_est, arma::mat alpha_est, arma::mat theta_est, const double gamma_est, arma::mat z_est, arma::mat w_est,
                              const double missing){
  double log_likelihood = 0, dist_temp = 0;
  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  int i, j, k;
  arma::dmat dist(nsample,nitem,fill::zeros);

  dist.fill(0.0);
  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      dist_temp = 0.0;
      for(j = 0; j < ndim; j++) dist_temp += std::pow((z_est(k,j)-w_est(i,j)), 2.0);
      dist(k,i) = std::sqrt(dist_temp);
    }
  }


  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      if(data(k,i) != missing){
        if(data(k,i) == 1.0) log_likelihood += -std::log(1.0 + std::exp(-(beta_est(i) + alpha_est(i) * theta_est(k) - gamma_est * dist(k,i))));
        else log_likelihood += -std::log(1.0 + std::exp(beta_est(i) + alpha_est(i) * theta_est(k) - gamma_est * dist(k,i)));
      }
    }
  }


  Rcpp::List output;
  output["log_likelihood"] = log_likelihood;

  return(output);

} // function end
