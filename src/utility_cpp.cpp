// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>

using namespace std;
using namespace arma;

// [[Rcpp::export]]
double Kumulavsech(arma::mat CC, double omega, arma::vec xlim, arma::vec ylim){
  double S = 0;
  int n = CC.n_rows;
  arma::vec s(n, fill::zeros);

  for(int i = 0; i < n; i++){
    s(i) = (normcdf(xlim(1), CC(i,0), omega) - normcdf(xlim(0), CC(i,0), omega))
    * (normcdf(ylim(1), CC(i,1), omega) - normcdf(ylim(0), CC(i,1), omega));
  }
  S = sum(s);

  return(S);
}


// [[Rcpp::export]]
double logpXCbeta(arma::mat X, arma::mat CC, double alpha, double omega, double AreaW, double integral){
  int nXrows = X.n_rows;
  int nCCrows = CC.n_rows;
  rowvec d;
  arma::vec res(nXrows, fill::zeros);
  arma::vec dummy;

  for(int i = 0; i < nXrows; i++){
    for(int j = 0; j < nCCrows; j++){
      d =  X.row(i) - CC.row(j);
      dummy =  (d*trans(d));
      res[i] = res[i] + exp(-dummy[0]/(2*omega*omega));
    }
  }
  double lik = AreaW - alpha*integral + nXrows*log(alpha/(2*(M_PI)*omega*omega)) + sum(log(res));

  return(lik);
}

// [[Rcpp::export]]
arma::vec NewPoint(arma::vec xlim, arma::vec ylim){
  double x, y = 0.0;
  arma::vec P(size(xlim), fill::zeros);

  x = randu(), y = randu();
  x = (xlim(1) - xlim(0)) * x + xlim(0);
  y = (ylim(1) - ylim(0)) * y + ylim(0);
  P(0) = x, P(1) = y;

  return(P);
}

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

static inline double inv_logit_stable_grm_ll(const double x){
  if(x >= 0.0){
    return 1.0 / (1.0 + std::exp(-x));
  }
  const double ex = std::exp(x);
  return ex / (1.0 + ex);
}

static inline double log_grm_prob_y0based_ll(const int y, const arma::rowvec& beta_i, const double eta){
  const int Kminus1 = (int)beta_i.n_elem;
  const int K = Kminus1 + 1;
  const double eps = 1e-16;

  auto p_ge = [&](const int k)->double{
    return inv_logit_stable_grm_ll(eta + beta_i((unsigned)(k-1)));
  };

  double p = 0.0;
  if(y <= 0){
    p = 1.0 - p_ge(1);
  }else if(y >= (K-1)){
    p = p_ge(K-1);
  }else{
    p = p_ge(y) - p_ge(y+1);
  }

  if(!std::isfinite(p) || p <= eps) p = eps;
  if(p >= 1.0 - eps) p = 1.0 - eps;
  return std::log(p);
}

// [[Rcpp::export]]
Rcpp::List log_likelihood_grm_cpp(arma::mat data, const int ndim, const int ncat,
                                  arma::mat beta_est, arma::mat theta_est,
                                  const double gamma_est, arma::mat z_est, arma::mat w_est,
                                  const double missing){
  if(ncat < 2){
    Rcpp::stop("ncat must be >= 2.");
  }

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  const int K = ncat;
  const int Kminus1 = K - 1;
  double log_likelihood = 0.0, dist_temp = 0.0;
  int i, j, k;

  // detect if data is 1-based (Likert 1..K) or 0-based (0..K-1)
  const double data_min = data.min();
  const int data_is_one_based = (data_min >= 1.0);

  if(beta_est.n_rows != nitem || beta_est.n_cols != Kminus1){
    Rcpp::stop("beta_est must be nitem x (ncat-1).");
  }

  arma::dmat dist(nsample, nitem, fill::zeros);
  dist.fill(0.0);
  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      dist_temp = 0.0;
      for(j = 0; j < ndim; j++) dist_temp += std::pow((z_est(k,j) - w_est(i,j)), 2.0);
      dist(k,i) = std::sqrt(dist_temp);
    }
  }

  for(i = 0; i < nitem; i++){
    arma::rowvec beta_i = beta_est.row(i);
    for(k = 0; k < nsample; k++){
      if(data(k,i) == missing) continue;
      int y = (int)(data(k,i));
      if(data_is_one_based) y -= 1;
      if(y < 0) y = 0;
      if(y > (K-1)) y = K-1;
      const double eta = theta_est(k) - gamma_est * dist(k,i);
      log_likelihood += log_grm_prob_y0based_ll(y, beta_i, eta);
    }
  }

  Rcpp::List output;
  output["log_likelihood"] = log_likelihood;
  return output;
}

// [[Rcpp::export]]
Rcpp::List log_likelihood_grm2pl_cpp(arma::mat data, const int ndim, const int ncat,
                                     arma::mat beta_est, arma::mat alpha_est, arma::mat theta_est,
                                     const double gamma_est, arma::mat z_est, arma::mat w_est,
                                     const double missing){
  if(ncat < 2){
    Rcpp::stop("ncat must be >= 2.");
  }

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  const int K = ncat;
  const int Kminus1 = K - 1;
  double log_likelihood = 0.0, dist_temp = 0.0;
  int i, j, k;

  // detect if data is 1-based (Likert 1..K) or 0-based (0..K-1)
  const double data_min = data.min();
  const int data_is_one_based = (data_min >= 1.0);

  if(beta_est.n_rows != nitem || beta_est.n_cols != Kminus1){
    Rcpp::stop("beta_est must be nitem x (ncat-1).");
  }
  if(alpha_est.n_rows != nitem){
    Rcpp::stop("alpha_est must have length nitem.");
  }

  arma::dmat dist(nsample, nitem, fill::zeros);
  dist.fill(0.0);
  for(i = 0; i < nitem; i++){
    for(k = 0; k < nsample; k++){
      dist_temp = 0.0;
      for(j = 0; j < ndim; j++) dist_temp += std::pow((z_est(k,j) - w_est(i,j)), 2.0);
      dist(k,i) = std::sqrt(dist_temp);
    }
  }

  for(i = 0; i < nitem; i++){
    arma::rowvec beta_i = beta_est.row(i);
    const double ai = alpha_est(i);
    for(k = 0; k < nsample; k++){
      if(data(k,i) == missing) continue;
      int y = (int)(data(k,i));
      if(data_is_one_based) y -= 1;
      if(y < 0) y = 0;
      if(y > (K-1)) y = K-1;
      const double eta = ai * theta_est(k) - gamma_est * dist(k,i);
      log_likelihood += log_grm_prob_y0based_ll(y, beta_i, eta);
    }
  }

  Rcpp::List output;
  output["log_likelihood"] = log_likelihood;
  return output;
}

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

#include <iostream>
#include <cstdio>
#include <chrono>
#include <Rcpp.h>

using namespace Rcpp;

void progressbar(int step, int total)
{
  // progress width
  const int pwidth = 72;
  
  // Prevent division by zero or incorrect display at the very beginning
  if (step <= 0 || total <= 0) return;
  if (step > total) step = total; // Ensure step doesn't exceed total
  
  // minus label len
  int pos = (static_cast<long long>(step) * pwidth) / total; // Use long long for intermediate calc
  int percent = (static_cast<long long>(step) * 100) / total;
  
  // calculate elapsed time
  auto current_time = std::chrono::system_clock::now();
  static auto start_time = current_time; // static to keep the initial value across calls *within the same load*
  if (step == 1) { // Reset start time when starting a new progress sequence
    start_time = current_time;
  }
  auto elapsed_duration = current_time - start_time;
  // Use milliseconds for potentially more precision in ETA for short tasks
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_duration).count();
  
  // calculate estimated time of arrival (ETA)
  long long remaining_ms = 0;
  if (step > 0 && elapsed_ms > 0) {
    // Estimate remaining time based on average time per step so far
    remaining_ms = (static_cast<long long>(total - step) * elapsed_ms) / step;
  }
  long long remaining_seconds_total = remaining_ms / 1000;
  
  // Calculate remaining time in hours, minutes, and seconds
  int hours = remaining_seconds_total / 3600;
  int minutes = (remaining_seconds_total % 3600) / 60;
  int seconds = remaining_seconds_total % 60;
  
  // Use Rprintf for output compatible with R console/RStudio
  Rprintf("[");
  for (int i = 0; i < pos; ++i) Rprintf("=");
  Rprintf("%*s", pwidth - pos, ""); // Print spaces to fill the bar
  Rprintf("] %3d%% ", percent);
  
  // Only print ETA if it's meaningful (elapsed time > 0)
  if (elapsed_ms > 0 && step < total) {
    Rprintf("ETA: %02d:%02d:%02d", hours, minutes, seconds);
  } else if (step == total) {
    // Optionally print total time instead of ETA when done
    long long elapsed_seconds_total = elapsed_ms / 1000;
    int elapsed_h = elapsed_seconds_total / 3600;
    int elapsed_m = (elapsed_seconds_total % 3600) / 60;
    int elapsed_s = elapsed_seconds_total % 60;
    Rprintf("Total: %02d:%02d:%02d", elapsed_h, elapsed_m, elapsed_s);
  } else {
    Rprintf("ETA: --:--:--"); // Indicate unknown ETA at the beginning
  }
  
  // Use \r (carriage return) to move cursor to the beginning of the line
  // This overwrites the previous progress bar line.
  Rprintf("\r");
  
  // Note: R's console often handles flushing automatically with \r or \n.
  // Explicit flushing (like fflush(stdout)) is usually not needed with Rprintf.
}
