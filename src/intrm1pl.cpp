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
Rcpp::List intrm1pl_cpp(arma::mat data, 
                  const int niter, const int nburn, const int nthin, const int nprint,
                  const double jump_beta, const double jump_theta, const double jump_delta,
                  const double pr_mean_beta, const double pr_sd_beta, const double pr_mean_theta,
                  const double pr_mean_delta, const double pr_sd_delta, const double pr_a_theta, const double pr_b_theta){

  int i, j, k, count, accept;
  double old_like_beta, new_like_beta, old_like_theta, new_like_theta, old_like_delta, new_like_delta, pr_sd_theta = 1.0;
  double num, den, ratio, un, post_a, post_b;
  
  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  
  arma::dvec oldbeta(nitem, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  arma::dvec newbeta = oldbeta;
  
  arma::dvec oldtheta(nsample, fill::randu);
  oldtheta = oldtheta * 4.0 - 2.0;
  arma::dvec newtheta = oldtheta;
  
  arma::dmat olddelta(nsample,nitem,fill::randu);
  olddelta = olddelta * 2.0 - 1.0;
  arma::dmat newdelta = olddelta;
  
  arma::dmat samp_beta((niter-nburn)/nthin, nitem, fill::zeros);
  arma::dmat samp_theta((niter-nburn)/nthin, nsample, fill::zeros);
  arma::dcube samp_delta(((niter-nburn)/nthin), nsample, nitem, fill::zeros);
  arma::dvec sample_sd_theta((niter-nburn)/nthin, fill::zeros);
  arma::dvec sample_mle((niter-nburn)/nthin, fill::zeros);
  arma::dvec accept_beta(nitem, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  arma::dmat accept_delta(nsample, nitem, fill::zeros);

  accept = count = 0;

  for(int iter = 0; iter < niter; iter++){

    // beta update
    for(i = 0; i < nitem; i++){
      newbeta(i) = R::rnorm(oldbeta(i), jump_beta);
      old_like_beta = new_like_beta = 0.0;
      for(k = 0; k < nsample; k++){
        if(data(k,i) == 1.0) old_like_beta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldtheta(k) + olddelta(k,i))));
        else old_like_beta += -std::log(1.0 + std::exp(oldbeta(i) + oldtheta(k) + olddelta(k,i)));
        if(data(k,i) == 1.0) new_like_beta += -std::log(1.0 + std::exp(-(newbeta(i) + oldtheta(k) + olddelta(k,i))));
        else new_like_beta += -std::log(1.0 + std::exp(newbeta(i) + oldtheta(k) + olddelta(k,i)));
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

    // theta update
    for(k = 0; k < nsample; k++){
      newtheta(k) = R::rnorm( oldtheta(k), jump_theta);
      old_like_theta = new_like_theta = 0.0;

      for(i = 0; i < nitem; i++){
        if(data(k,i) == 1.0) old_like_theta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldtheta(k) + olddelta(k,i))));
        else old_like_theta += -std::log(1.0 + std::exp(oldbeta(i) + oldtheta(k)+ olddelta(k,i)));
        if(data(k,i) == 1.0) new_like_theta += -std::log(1.0 + std::exp(-(oldbeta(i) + newtheta(k) + olddelta(k,i))));
        else new_like_theta += -std::log(1.0 + std::exp(oldbeta(i) + newtheta(k)+ olddelta(k,i)));
      }
      num = new_like_theta + R::dnorm4(newtheta(k), pr_mean_theta, pr_sd_theta, 1);
      den = old_like_theta + R::dnorm4(oldtheta(k), pr_mean_theta, pr_sd_theta, 1);
      ratio = num - den;

      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }

      if(accept == 1){
        oldtheta(k) = newtheta(k);
        accept_theta(k) += 1.0 / (niter * 1.0);
      }
      else newtheta(k) = oldtheta(k);
    }

    // dela update
    for(k = 0; k < nsample; k++){
      for(i = 0; i < nitem; i++){
        newdelta(k,i) = R::rnorm(olddelta(k,i), jump_delta);
        old_like_delta = new_like_delta = 0.0;
        
        if(data(k,i) == 1.0) new_like_delta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldtheta(k) + newdelta(k,i))));
        else new_like_delta += -std::log(1.0 + std::exp(oldbeta(i) + oldtheta(k) + newdelta(k,i)));
        if(data(k,i) == 1.0) old_like_delta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldtheta(k) + olddelta(k,i))));
        else old_like_delta += -std::log(1.0 + std::exp(oldbeta(i) + oldtheta(k) + olddelta(k,i)));
        
        num = new_like_delta + R::dnorm4(newdelta(k,i),pr_mean_delta,pr_sd_delta,1);
        den = old_like_delta + R::dnorm4(olddelta(k,i),pr_mean_delta,pr_sd_delta,1);
        ratio = num - den;
        
        if(ratio > 0.0) accept = 1;
        else{
          un = R::runif(0,1);
          if(std::log(un) < ratio) accept = 1;
          else accept = 0;
        }
        
        if(accept == 1){
          olddelta(k,i) = newdelta(k,i);
          accept_delta(k,i) += 1.0 / (niter * 1.0);
        }
        else{
          newdelta(k,i) = olddelta(k,i);
        }
      }
    }
    

    //sigma_theta update with gibbs
    post_a = 2 * pr_a_theta  + nsample;
    post_b = pr_b_theta;
    for(j = 0; j < nsample; j++) post_b += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2;
    pr_sd_theta = std::sqrt(2 * post_b *(1.0 /  R::rchisq(post_a)));

    
    //Burn in 
    if(iter >= nburn && iter % nthin == 0){
      for(i = 0; i < nitem; i++) samp_beta(count,i) = oldbeta(i);
      for(k = 0; k < nsample; k++) samp_theta(count,k) = oldtheta(k);
      for(k = 0; k < nsample; k++)
        for(i = 0; i < nitem; i++) samp_delta(count,k,i) = olddelta(k,i);
      sample_sd_theta(count) = pr_sd_theta;
      count++;
    }

    if(iter % nprint == 0){
      Rprintf("Iteration: %.5u", iter); 
      for(i = 0 ; i < nitem ; i++ ) {
        Rprintf("% .3f ", oldbeta(i));
      }
      Rprintf("%.3f\n", pr_sd_theta);
    }
  }
  
  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["delta"] = samp_delta;
  output["sigma_theta"] = sample_sd_theta;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  output["accept_delta"] = accept_delta;
  
  return(output);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

