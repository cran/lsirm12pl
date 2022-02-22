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
Rcpp::List intrm2pl_mar_cpp(arma::mat data, const int niter, const int nburn, const int nthin, const int nprint,
                        const double jump_beta, const double jump_theta, const double jump_alpha, const double jump_delta,
                        const double pr_mean_beta, const double pr_sd_beta, const double pr_mean_theta,
                        const double pr_mean_delta, const double pr_sd_delta, const double pr_mean_alpha, const double pr_sd_alpha, 
                        const double pr_a_theta, const double pr_b_theta, const double missing){

  int i, j, k, count, accept;
  double old_like_beta, new_like_beta, old_like_theta, new_like_theta, old_like_delta, new_like_delta, pr_sd_theta = 1.0,
    old_like_alpha, new_like_alpha ;
  double num, den, ratio, un, post_a, post_b;

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  
  arma::dvec oldbeta(nitem, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  arma::dvec newbeta = oldbeta;
  
  arma::dvec oldtheta(nsample, fill::randu);
  oldtheta = oldtheta * 4.0 - 2.0;
  arma::dvec newtheta = oldtheta;
  
  arma::dvec oldalpha(nitem , fill::randu);
  // oldalpha = oldalpha * 4.0 - 2.0;
  oldalpha(0) = 1;
  arma::dvec newalpha = oldalpha;
  
  arma::dmat olddelta(nsample,nitem,fill::randu);
  olddelta = olddelta * 2.0 - 1.0;
  arma::dmat newdelta = olddelta;
  
  arma::dmat samp_beta((niter-nburn)/nthin-1, nitem, fill::zeros);
  arma::dmat samp_theta((niter-nburn)/nthin-1, nsample, fill::zeros);
  arma::dmat samp_alpha((niter-nburn)/nthin-1, nitem, fill::zeros);
  arma::dcube samp_delta(((niter-nburn)/nthin), nsample, nitem, fill::zeros);
  arma::dvec sample_sd_theta((niter-nburn)/nthin-1, fill::zeros);
  arma::dvec accept_beta(nitem, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  arma::dvec accept_alpha(nitem, fill::zeros);
  arma::dmat accept_delta(nsample, nitem, fill::zeros);
  
  accept = count = 0;

  double p_ki, impute_value;
  int nmissing, mi;
  //missing_data is indicator matrix for missing
  arma::mat missing_data(data.begin(), data.n_rows, data.n_cols, true);
  nmissing = 0;
  for(i =0; i<nitem; i++){
    for(k =0; k<nsample; k++){
      if(missing_data(k,i) == missing){
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
      for(k =0; k<nsample; k++){
        if(missing_data(k,i) == missing){
          p_ki = 1 / (1 + std::exp(-(oldbeta(i) + oldalpha(i) * oldtheta(k) + olddelta(k,i))));
          impute_value = R::rbinom(1, p_ki);
          data(k,i) = impute_value;
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
        if(data(j,i) == 1.0) new_like_beta += -std::log(1.0 + std::exp(-(newbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i)))) ;
        else new_like_beta += -std::log(1.0 + std::exp(newbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i))) ;
        if(data(j,i) == 1.0) old_like_beta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i)))) ;
        else old_like_beta += -std::log(1.0 + std::exp(oldbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i))) ;
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
        if(data(j,i) == 1.0) new_like_theta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldalpha(i) * newtheta(j) + olddelta(j,i)))) ;
        else new_like_theta += -std::log(1.0 + std::exp(oldbeta(i) + oldalpha(i) * newtheta(j) + olddelta(j,i))) ;
        if(data(j,i) == 1.0) old_like_theta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i)))) ;
        else old_like_theta += -std::log(1.0 + std::exp(oldbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i))) ;
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
    
    // alpha update
    for(i = 1; i < nitem; i++){
      newalpha(i) = R::rlnorm(std::log(oldalpha(i)) , jump_alpha);
      old_like_alpha = new_like_alpha = 0.0;
      
      for(k = 0; k < nsample; k++){
        if(data(k,i) == 1.0) new_like_beta += -std::log(1.0 + std::exp(-(oldbeta(i) + newalpha(i) * oldtheta(k) + olddelta(k,i))));
        else new_like_beta += -std::log(1.0 + std::exp(oldbeta(i) + newalpha(i) * oldtheta(k) + olddelta(k,i)));
        if(data(k,i) == 1.0) old_like_beta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldalpha(i) * oldtheta(k) + olddelta(k,i))));
        else old_like_beta += -std::log(1.0 + std::exp(oldbeta(i) + oldalpha(i) * oldtheta(k) + olddelta(k,i)));
      }
      
      num = new_like_alpha + R::dlnorm(oldalpha(i), std::log(newalpha(i)) , jump_alpha, 1) + R::dlnorm(newalpha(i), pr_mean_alpha, pr_sd_alpha, 1);
      den = old_like_alpha + R::dlnorm(newalpha(i), std::log(oldalpha(i)) , jump_alpha, 1) + R::dlnorm(oldalpha(i), pr_mean_alpha, pr_sd_alpha, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldalpha(i) = newalpha(i);
        accept_alpha(i) += 1.0 / (niter * 1.0);
      }
      else newalpha(i) = oldalpha(i);    
    }
    
    // delta update
    for(j = 0; j < nsample; j++){
      for(i = 0; i < nitem; i++){
        newdelta(j,i) = R::rnorm(olddelta(j,i), jump_delta);
        old_like_delta = new_like_delta = 0.0;
        
        if(data(j,i) == 1.0) new_like_delta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldalpha(i) * oldtheta(j) + newdelta(j,i))));
        else new_like_delta += -std::log(1.0 + std::exp(oldbeta(i) + oldalpha(i) * oldtheta(j) + newdelta(j,i)));
        if(data(j,i) == 1.0) old_like_delta += -std::log(1.0 + std::exp(-(oldbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i))));
        else old_like_delta += -std::log(1.0 + std::exp(oldbeta(i) + oldalpha(i) * oldtheta(j) + olddelta(j,i)));
        
        num = new_like_delta + R::dnorm4(newdelta(j,i),pr_mean_delta,pr_sd_delta,1);
        den = old_like_delta + R::dnorm4(olddelta(j,i),pr_mean_delta,pr_sd_delta,1);
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
    post_a = 2 * pr_a_theta  + nsample;
    post_b = pr_b_theta;
    for(j = 0; j < nsample; j++) post_b += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2;
    pr_sd_theta = std::sqrt(2 * post_b *(1.0 /  R::rchisq(post_a)));
    
    //burn in
    if(iter > nburn && iter % nthin == 0){
      for(i = 0; i < nitem; i++) samp_beta(count,i) = oldbeta(i);
      for(i = 0; i < nitem; i++) samp_alpha(count,i) = oldalpha(i);
      for(j = 0; j < nsample; j++) samp_theta(count,j) = oldtheta(j);
      for(j = 0; j < nsample; j++)
        for(i = 0; i < nitem; i++) samp_delta(count,j,i) = olddelta(j,i);
      sample_sd_theta(count) = pr_sd_theta;
      for(mi = 0; mi < nmissing; mi++) samp_impute(count,mi) = impute_col(mi);          
      count++;
    }
    
    if(iter % nprint == 0){
      Rprintf("Iteration: %.5u ", iter); 
      for(i = 0 ; i < nitem ; i++ ) {
        Rprintf("% .3f ", oldbeta(i));
      }
      Rprintf(" %.3f\n", pr_sd_theta);
    }
  }
  
  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["alpha"] = samp_alpha;
  output["delta"] = samp_delta;
  output["sigma_theta"] = sample_sd_theta;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  output["accept_alpha"] = accept_alpha;
  output["accept_delta"] = accept_delta;
  output["impute"] = samp_impute;    
  
  return(output);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

