// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
// [[Rcpp::export]]
Rcpp::List onepl_cpp(arma::mat data, const int niter, const int nburn, const int nthin, const int nprint,
                  const double jump_beta, const double jump_theta, 
                  const double pr_mean_beta, const double pr_sd_beta, const double pr_mean_theta, 
                  const double pr_a_theta, const double pr_b_theta){

  int i, j, count, accept;
  double old_like_beta, new_like_beta, old_like_theta, new_like_theta, pr_sd_theta = 1.0;
  double num, den, ratio, un, post_a, post_b;

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  
  arma::dvec oldbeta(nitem, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  arma::dvec newbeta = oldbeta;
  
  arma::dvec oldtheta(nsample, fill::randu);
  oldtheta = oldtheta * 4.0 - 2.0;
  arma::dvec newtheta = oldtheta;
  
  arma::dvec samp_like(nsample, fill::zeros);
  arma::dvec item_like(nsample, fill::zeros);
  
  arma::dmat samp_beta((niter-nburn)/nthin-1, nitem, fill::zeros);
  arma::dmat samp_theta((niter-nburn)/nthin-1, nsample, fill::zeros);
  arma::dvec sample_sd_theta((niter-nburn)/nthin-1, fill::zeros);
  arma::dvec accept_beta(nitem, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  
  accept = count = 0;
  for(int iter = 0; iter < niter; iter++){
    for(i = 0; i < nitem; i++){
      newbeta(i) = R::rnorm( oldbeta(i), jump_beta);
      old_like_beta = new_like_beta = 0.0;
      
      for(int k = 0; k < nsample; k++){
        if(data(k,i) == 1.0) old_like_beta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldtheta(k))));
        else old_like_beta += -std::log(1.0 + std::exp(oldbeta(i) + oldtheta(k)));
        if(data(k,i) == 1.0) new_like_beta += -std::log(1.0 + std::exp(-(newbeta(i) + oldtheta(k))));
        else new_like_beta += -std::log(1.0 + std::exp(newbeta(i) + oldtheta(k)));
      }
      
      num = new_like_beta + R::dnorm4(newbeta(i), pr_mean_beta, pr_sd_beta, 1);
      den = old_like_beta + R::dnorm4(oldbeta(i), pr_mean_beta, pr_sd_beta, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldbeta(i) = newbeta(i);
        accept_beta(i) += 1.0 / (niter * 1.0);
      }
      else newbeta(i) = oldbeta(i);
      
    }
    
    for(j = 0; j < nsample; j++){
      newtheta(j) = R::rnorm(oldtheta(j), jump_theta);
      old_like_theta = new_like_theta = 0.0;
      
      for(int k = 0; k < nitem; k++){
        if(data(j,k) == 1.0) old_like_theta += -std::log(1.0 + std::exp(-(oldbeta(k) + oldtheta(j))));
        else old_like_theta += -std::log(1.0 + std::exp(oldbeta(k) + oldtheta(j)));
        if(data(j,k) == 1.0) new_like_theta += -std::log(1.0 + std::exp(-(oldbeta(k) + newtheta(j))));
        else new_like_theta += -std::log(1.0 + std::exp(oldbeta(k) + newtheta(j)));
      }
      num = new_like_theta + R::dnorm4(newtheta(j), pr_mean_theta, pr_sd_theta, 1);
      den = old_like_theta + R::dnorm4(oldtheta(j), pr_mean_theta, pr_sd_theta, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldtheta(j) = newtheta(j);
        accept_theta(j) += 1.0 / (niter * 1.0);
      }
      else newtheta(j) = oldtheta(j);
    }
    
    //sigma_theta update with gibbs
    post_a = 2 * pr_a_theta  + nsample;
    post_b = pr_b_theta;
    for(j = 0; j < nsample; j++) post_b += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2;
    pr_sd_theta = std::sqrt(2 * post_b *(1.0 /  R::rchisq(post_a)));
    
    if(iter > nburn && iter % nthin == 0){
      for(int i = 0; i < nitem; i++) samp_beta(count,i) = oldbeta(i);
      for(int k = 0; k < nsample; k++) samp_theta(count, k) = oldtheta(k);
      sample_sd_theta(count) = pr_sd_theta;
      count++;
    }
  
    if(iter % nprint == 0){
      Rprintf("Iteration: %.5u", iter); 
      for(i = 0 ; i < 6 ; i++ ) {
        Rprintf("% .3f ", oldbeta(i));
      }
      Rprintf("%.3f\n", pr_sd_theta);
    }
  }
  
  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["sigma_theta"] = sample_sd_theta;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  
  return(output);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

