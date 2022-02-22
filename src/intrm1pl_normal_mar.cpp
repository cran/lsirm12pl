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
Rcpp::List intrm1pl_normal_mar_cpp(arma::mat data, const int niter, const int nburn, const int nthin, const int nprint,
                             const double jump_beta, const double jump_theta, const double jump_delta,
                             const double pr_mean_beta, const double pr_sd_beta, const double pr_mean_theta, const double pr_mean_delta, const double pr_sd_delta,
                             const double pr_a_eps, const double pr_b_eps, const double pr_a_theta, const double pr_b_theta, const double missing){

  int i, j, count, accept;
  double old_like_beta, new_like_beta, old_like_theta, new_like_theta, old_like_delta, new_like_delta, pr_sd_theta = 1.0, pr_sd ;
  double num, den, ratio, un, post_a_sigma, post_b_sigma, post_a_th_sigma, post_b_th_sigma;
  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  pr_sd = 1;
  
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
  arma::dvec sample_sd((niter-nburn)/nthin, fill::zeros);
  arma::dvec accept_beta(nitem, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  arma::dmat accept_delta(nsample, nitem, fill::zeros);
  
  accept = count = 0;

  double impute_mean, impute_value;
  int nmissing, mi;
  //missing_data is indicator matrix for missing
  arma::mat missing_data(data.begin(), data.n_rows, data.n_cols, true);
  nmissing = 0;
  for(i =0; i<nitem; i++){
    for(j =0; j<nsample; j++){
      if(missing_data(j,i) == missing){
        nmissing++;
      }
    }
  }
  arma::dvec impute_col(nmissing, fill::randu);
  arma::dmat samp_impute((niter-nburn)/nthin, nmissing, fill::zeros); 

  for(int iter = 0; iter < niter; iter++){
    // Imputation step
    mi = 0;
    for(i =0; i<nitem; i++){
      for(j =0; j<nsample; j++){
        if(missing_data(j,i) == missing){
          impute_mean = oldbeta(i) + oldtheta(j) + olddelta(j,i);
          impute_value = R::rnorm( impute_mean, pr_sd);
          data(j,i) = impute_value;
          impute_col(mi) = impute_value;
          mi++;
        }
      }
    }                      
    
    //beta update
    for(i = 0; i < nitem; i++){
      newbeta(i) = R::rnorm( oldbeta(i), jump_beta);
      old_like_beta = new_like_beta = 0.0;
     
     for(j = 0; j < nsample; j++){
       new_like_beta += - pow((data(j,i) - newbeta(i) - oldtheta(j) - olddelta(j,i)), 2) / (2 * pow(pr_sd, 2)) ;
       old_like_beta += - pow((data(j,i) - oldbeta(i) - oldtheta(j) - olddelta(j,i)), 2) / (2 * pow(pr_sd, 2)) ;
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
    for(j = 0; j < nsample; j++){
      newtheta(j) = R::rnorm(oldtheta(j), jump_theta);
      old_like_theta = new_like_theta = 0.0;
      
      for(i = 0; i < nitem; i++){
        new_like_theta += - pow((data(j,i) - oldbeta(i) - newtheta(j) - olddelta(j,i)),2) / (2 * pow(pr_sd, 2)) ;
        old_like_theta += - pow((data(j,i) - oldbeta(i) - oldtheta(j) - olddelta(j,i)),2) / (2 * pow(pr_sd, 2)) ;
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
    
    // delta update
    for(j = 0; j < nsample; j++){
      for(i = 0; i < nitem; i++){
        newdelta(j,i) = R::rnorm(olddelta(j,i), jump_delta);
        old_like_delta = new_like_delta = 0.0;
        
        new_like_delta += - pow((data(j,i) - oldbeta(i) - oldtheta(j) - newdelta(j,i)), 2) / (2 * pow(pr_sd, 2)) ;
        old_like_delta += - pow((data(j,i) - oldbeta(i) - oldtheta(j) - olddelta(j,i)), 2)  / (2 * pow(pr_sd, 2)) ;
        
        num = new_like_delta + R::dnorm4(newdelta(j,i), pr_mean_delta, pr_sd_delta, 1);
        den = old_like_delta + R::dnorm4(olddelta(j,i), pr_mean_delta, pr_sd_delta, 1);
        ratio = num - den;
        
        if(ratio > 0.0) accept = 1;
        else{
          un = R::runif(0,1);
          if(std::log(un) < ratio) accept = 1;
          else accept = 0;
        }
        
        if(accept == 1){
          olddelta(j,i) = newdelta(j,i);
          accept_delta(j,i) += 1.0 / (niter * 1.0);
        }
        else{
          newdelta(j,i) = olddelta(j,i);
        }
      }
    }
    
    
    //sigma_theta update with gibbs
    post_a_th_sigma = 2 * pr_a_theta  + nsample;
    post_b_th_sigma = pr_b_theta;
    for(j = 0; j < nsample; j++) post_b_th_sigma += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2;
    pr_sd_theta = std::sqrt(2 * post_b_th_sigma *(1.0 /  R::rchisq(post_a_th_sigma)));
    
    //sigma update with gibbs
    
    post_a_sigma = 2 * pr_a_eps  + nsample * nitem ;
    post_b_sigma = pr_b_eps;
    for(j = 0; j < nsample; j++){
      for(i = 0; i < nitem; i++) post_b_sigma += std::pow((data(j,i) - oldbeta(i) - oldtheta(j) - olddelta(j,i)), 2.0) / 2;
    }
    pr_sd = std::sqrt(2 * post_b_sigma *(1.0 /  R::rchisq(post_a_sigma)));
    
    //burn in
    if(iter >= nburn && iter % nthin == 0){
      for(i = 0; i < nitem; i++) samp_beta(count,i) = oldbeta(i);
      for(j = 0; j < nsample; j++) samp_theta(count,j) = oldtheta(j);
      for(j = 0; j < nsample; j++)
        for(i = 0; i < nitem; i++) samp_delta(count,j,i) = olddelta(j,i);
      sample_sd_theta(count) = pr_sd_theta;
      sample_sd(count) = pr_sd;
      for(mi = 0; mi < nmissing; mi++) samp_impute(count,mi) = impute_col(mi);        
      count++;
    }
    
    if(iter % nprint == 0){
      Rprintf("Iteration: %.5u ", iter); 
      for(i = 0 ; i < nitem ; i++ ) {
        Rprintf("% .3f ", oldbeta(i));
      }
      Rprintf(" %.3f\n", pr_sd_theta);
      Rprintf(" %.3f\n", pr_sd);
    }
  }
  
  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["delta"] = samp_delta;
  output["sigma_theta"] = sample_sd_theta;
  output["sigma"] = sample_sd;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  output["accept_delta"] = accept_delta;
  output["impute"] = samp_impute;    
  
  return(output);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

