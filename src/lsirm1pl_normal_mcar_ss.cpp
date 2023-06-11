// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include "progress.h"
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
Rcpp::List lsirm1pl_normal_mcar_ss_cpp(arma::mat data, const int ndim, const int niter, const int nburn, const int nthin, const int nprint,
                              const double jump_beta, const double jump_theta, const double jump_gamma, const double jump_z, const double jump_w,
                              const double pr_mean_beta, const double pr_sd_beta, const double pr_mean_theta,
                              const double pr_spike_mean, const double pr_spike_sd, const double pr_slab_mean, const double pr_slab_sd,
                              const double pr_a_theta, const double pr_b_theta, const double pr_a_eps, const double pr_b_eps,
                              const double pr_beta_a, const double pr_beta_b, const double missing, const bool verbose){

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;

  int i, j, k, count, accept, num_notmissing;
  double num, den, old_like_beta, new_like_beta, old_like_theta, new_like_theta;
  double old_like_z, new_like_z, old_like_w, new_like_w, old_like_gamma, new_like_gamma;
  double ratio, un, dist_temp, dist_old_temp, dist_new_temp;
  double post_a_sigma, post_b_sigma, post_a_th_sigma, post_b_th_sigma;
  double pr_mean_z = 0.0, pr_sd_z = 1.0, pr_mean_w = 0.0, pr_sd_w = 1.0, pr_sd = 1.0, pr_sd_theta = 1.0, mle;


  arma::dvec oldbeta(nitem, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  arma::dvec newbeta = oldbeta;

  arma::dvec oldtheta(nsample, fill::randu);
  oldtheta = oldtheta * 4.0 - 2.0;
  arma::dvec newtheta = oldtheta;

  arma::dmat oldz(nsample,ndim,fill::randu);
  oldz = oldz * 2.0 - 1.0;
  arma::dmat newz = oldz;

  arma::dmat oldw(nitem,ndim,fill::randu);
  oldw = oldw * 2.0 - 1.0;
  arma::dmat neww = oldw;


  double oldgamma = 1, newgamma = 1 , pi_prob = 0.5, xi;
  int pi;
  pi = R::rbinom(1, pi_prob);
  xi = R::rbeta(pr_beta_a + pi, pr_beta_b + 1.0 - pi);

  arma::dmat samp_beta((niter-nburn)/nthin, nitem, fill::zeros);
  arma::dmat samp_theta((niter-nburn)/nthin, nsample, fill::zeros);
  arma::dcube samp_z(((niter-nburn)/nthin), nsample, ndim, fill::zeros);
  arma::dcube samp_w(((niter-nburn)/nthin), nitem, ndim, fill::zeros);
  arma::dvec samp_sd_theta((niter-nburn)/nthin, fill::zeros);
  arma::dvec samp_sd((niter-nburn)/nthin, fill::zeros);
  arma::dvec samp_mle((niter-nburn)/nthin, fill::zeros);
  arma::dvec samp_gamma(((niter-nburn)/nthin), fill::zeros);
  arma::ivec samp_pi(((niter-nburn)/nthin), fill::zeros);
  arma::dvec samp_xi(((niter-nburn)/nthin), fill::zeros);

  arma::dvec accept_beta(nitem, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  arma::dvec accept_z(nsample, fill::zeros);
  arma::dvec accept_w(nitem, fill::zeros);
  double accept_gamma=0;

  accept = count = 0;

  arma::dmat dist(nsample,nitem,fill::zeros);
  arma::dvec old_dist_k(nitem,fill::zeros);
  arma::dvec new_dist_k(nitem,fill::zeros);
  arma::dvec old_dist_i(nsample,fill::zeros);
  arma::dvec new_dist_i(nsample,fill::zeros);

  for(int iter = 0; iter < niter; iter++){
      if (iter % 10 == 0){
      Rcpp::checkUserInterrupt();
    }


    //dist(j,i) is distance of z_j and w_i
    dist.fill(0.0);
    for(i = 0; i < nitem; i++){
      for(k = 0; k < nsample; k++){
        dist_temp = 0.0;
        for(j = 0; j < ndim; j++) dist_temp += std::pow((oldz(k,j)-oldw(i,j)), 2.0);
        dist(k,i) = std::sqrt(dist_temp);
      }
    }

    // beta update
    for(i = 0; i < nitem; i++){
      newbeta(i) = R::rnorm(oldbeta(i), jump_beta);
      old_like_beta = new_like_beta = 0.0;
      for(k = 0; k < nsample; k++){
        if(data(k,i) != missing){
          new_like_beta += - pow((data(k,i) - newbeta(i) - oldtheta(k) + oldgamma * dist(k,i)), 2) / (2 * pow(pr_sd, 2)) ;
          old_like_beta += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * dist(k,i)), 2) / (2 * pow(pr_sd, 2)) ;
        }
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
        if(data(k,i) != missing){
        new_like_theta += - pow((data(k,i) - oldbeta(i) - newtheta(k) + oldgamma * dist(k,i)),2) / (2 * pow(pr_sd, 2)) ;
        old_like_theta += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * dist(k,i)),2) / (2 * pow(pr_sd, 2)) ;
        }
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

    // gamma(log(gamma)) update
    newgamma = R::rlnorm(std::log(oldgamma), jump_gamma);
    old_like_gamma = new_like_gamma = 0.0;

    for(k = 0; k < nsample; k++){
      for(i = 0; i < nitem; i++){
        if(data(k,i) != missing){
        new_like_gamma += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + newgamma * dist(k,i)),2) / (2 * pow(pr_sd, 2)) ;
        old_like_gamma += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * dist(k,i)),2) / (2 * pow(pr_sd, 2)) ;
        }
      }
    }

    num = new_like_gamma + R::dlnorm(oldgamma, std::log(newgamma), jump_gamma, 1) ;
    den = old_like_gamma + R::dlnorm(newgamma, std::log(oldgamma), jump_gamma, 1) ;
    if(pi == 1){
      num += R::dlnorm(newgamma, pr_slab_mean, pr_slab_sd, 1);
      den += R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 1);
    }
    else{
      num += R::dlnorm(newgamma, pr_spike_mean, pr_spike_sd, 1);
      den += R::dlnorm(oldgamma, pr_spike_mean, pr_spike_sd, 1);
    }

    ratio = num - den;

    if(ratio > 0.0) accept = 1;
    else{
      un = R::runif(0,1);
      if(std::log(un) < ratio) accept = 1;
      else accept = 0;
    }

    if(accept == 1){
      oldgamma = newgamma;
      accept_gamma += 1.0 / (niter * 1.0);
    }
    else newgamma = oldgamma;

    //pi update(spike slab parameter)
    pi_prob  = xi * R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 0);
    pi_prob /= (1.0 - xi) * R::dlnorm(oldgamma, pr_spike_mean, pr_spike_sd, 0) + xi * R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 0);
    pi = R::rbinom(1, pi_prob);

    // xi update(spike slab parameter)
    xi = R::rbeta(pr_beta_a + pi, pr_beta_b + 1.0 - pi);

    // zj update
    for(k = 0; k < nsample; k++){
      for(j = 0; j < ndim; j++) newz(k,j) = R::rnorm(oldz(k,j), jump_z);
      old_like_z = new_like_z = 0.0;

      //calculate distance of oldw and newz
      for(i = 0; i < nitem; i++){
        dist_old_temp = dist_new_temp = 0.0;
        for(j = 0; j < ndim; j++){
          dist_new_temp += std::pow((newz(k,j)-oldw(i,j)), 2.0);
          dist_old_temp += std::pow((oldz(k,j)-oldw(i,j)), 2.0);
        }
        new_dist_k(i) = sqrt(dist_new_temp);
        old_dist_k(i) = sqrt(dist_old_temp);
      }

      //calculate likelihood
      for(i = 0; i < nitem; i++){
        if(data(k,i) != missing){
        new_like_z += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * new_dist_k(i)),2) / (2 * pow(pr_sd, 2)) ;
        old_like_z += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * old_dist_k(i)),2) / (2 * pow(pr_sd, 2)) ;
        }
      }

      num = den = 0.0;
      for(j = 0; j < ndim; j++){
        num += R::dnorm4(newz(k,j),pr_mean_z,pr_sd_z,1);
        den += R::dnorm4(oldz(k,j),pr_mean_z,pr_sd_z,1);
      }
      //Rprintf("%.3f %.3f %.3f %.3f\n", num, den, new_like_z, old_like_z);
      //arma::dvec newzz = dmvnorm(newz.cols(2*j,2*j+1),pr_mean_z,pr_cov_z,TRUE);
      //arma::dvec oldzz = dmvnorm(oldz.cols(2*j,2*j+1),pr_mean_z,pr_cov_z,TRUE);

      num += new_like_z;
      den += old_like_z;
      ratio = num - den;

      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }

      if(accept == 1){
        for(j = 0; j < ndim; j++) oldz(k,j) = newz(k,j);
        accept_z(k) += 1.0 / (niter * 1.0);
      }
      else{
        for(j = 0; j < ndim; j++) newz(k,j) = oldz(k,j);
      }
    }

    // wi update
    for(i = 0; i < nitem; i++){
      for(j = 0; j < ndim; j++) neww(i,j) = R::rnorm(oldw(i,j), jump_w);
      old_like_w = new_like_w = 0.0;

      //calculate distance of neww and oldz
      for(k = 0; k < nsample; k++){
        dist_old_temp = dist_new_temp = 0.0;
        for(j = 0; j < ndim; j++){
          dist_new_temp += std::pow((oldz(k,j)-neww(i,j)), 2.0);
          dist_old_temp += std::pow((oldz(k,j)-oldw(i,j)), 2.0);
        }
        new_dist_i(k) = sqrt(dist_new_temp);
        old_dist_i(k) = sqrt(dist_old_temp);
      }

      //calculate likelihood
      for(k = 0; k < nsample; k++){
        if(data(k,i) != missing){
        new_like_w += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * new_dist_i(k)),2) / (2 * pow(pr_sd, 2)) ;
        old_like_w += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * old_dist_i(k)),2) / (2 * pow(pr_sd, 2)) ;
        }
      }

      num = den = 0.0;
      for(j = 0; j < ndim; j++){
        num += R::dnorm4(neww(i,j),pr_mean_w,pr_sd_w,1);
        den += R::dnorm4(oldw(i,j),pr_mean_w,pr_sd_w,1);
      }

      num += new_like_w;
      den += old_like_w;
      ratio = num - den;

      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }

      if(accept == 1){
        for(j = 0; j < ndim; j++) oldw(i,j) = neww(i,j);
        accept_w(i) += 1.0 / (niter * 1.0);
      }
      else{
        for(j = 0; j < ndim; j++) neww(i,j) = oldw(i,j);
      }
    }


    //sigma_theta update with gibbs
    post_a_th_sigma = 2 * pr_a_theta  + nsample;
    post_b_th_sigma = pr_b_theta;
    for(j = 0; j < nsample; j++) post_b_th_sigma += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2;
    pr_sd_theta = std::sqrt(2 * post_b_th_sigma *(1.0 /  R::rchisq(post_a_th_sigma)));

    //dist(j,i) is distance of z_j and w_i
    dist.fill(0.0);
    for(i = 0; i < nitem; i++){
      for(k = 0; k < nsample; k++){
        dist_temp = 0.0;
        for(j = 0; j < ndim; j++) dist_temp += std::pow((oldz(k,j)-oldw(i,j)), 2.0);
        dist(k,i) = std::sqrt(dist_temp);
      }
    }


    //sigma update with gibbs
    num_notmissing = 0;
    for(i = 0; i < nitem; i++){
      for(k = 0; k < nsample; k++){
        if(data(k,i) != missing){
          num_notmissing += 1;
        }
      }
    }
    post_a_sigma = 2 * pr_a_eps  + num_notmissing * nitem ;
    post_b_sigma = pr_b_eps;
    for(i = 0; i < nitem; i++){
      for(k = 0; k < nsample; k++){
        if(data(k, i) != missing){
          post_b_sigma += std::pow((data(k, i) - oldbeta(i) - oldtheta(k) + oldgamma * dist(k, i)), 2.0) / 2;
        }
      }
    }
    pr_sd = std::sqrt(2 * post_b_sigma *(1.0 /  R::rchisq(post_a_sigma)));

    // burn, thin
    if(iter >= nburn && iter % nthin == 0){
      for(i = 0; i < nitem; i++) samp_beta(count,i) = oldbeta(i);
      for(k = 0; k < nsample; k++) samp_theta(count,k) = oldtheta(k);
      for(i = 0; i < nitem; i++){
        for(j = 0; j < ndim; j++){
          samp_w(count,i,j) = oldw(i,j);
        }
      }
      for(k = 0; k < nsample; k++){
        for(j = 0; j < ndim; j++){
          samp_z(count,k,j) = oldz(k,j);
        }
      }

      samp_gamma(count) = oldgamma;
      samp_pi(count) = pi;
      samp_xi(count) = xi;
      samp_sd_theta(count) = pr_sd_theta;
      samp_sd(count) = pr_sd;

      mle = 0.0;
      for(i = 0; i < nitem; i++) mle += R::dnorm4(oldbeta(i), pr_mean_beta, pr_sd_beta, 1);
      for(k = 0; k < nsample; k++) mle += R::dnorm4(oldtheta(k), pr_mean_theta, pr_sd_theta, 1);
      for(i = 0; i < nitem; i++)
        for(j = 0; j < ndim; j++) mle += R::dnorm4(oldw(i,j),pr_mean_w,pr_sd_w,1);
      for(k = 0; k < nsample; k++)
        for(j = 0; j < ndim; j++) mle += R::dnorm4(oldz(k,j),pr_mean_z,pr_sd_z,1);
      for(k = 0; k < nsample; k++){
        for(i = 0; i < nitem; i++){
          if(data(k,i) != missing){
          mle += - pow((data(k,i) - oldbeta(i) - oldtheta(k) + oldgamma * dist(k,i)), 2) / (2 * pow(pr_sd, 2)) ;
          }
        }
      }
      if(pi == 1) mle += R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 1);
      else mle += R::dnorm(oldgamma, pr_spike_mean, pr_spike_sd, 1);

      samp_mle(count) = mle;

      count++;
    }
    
    if(verbose){
      int percent = 0;
      if(iter % nprint == 0){
        percent = (iter*100)/niter;
        Rprintf("Iteration: %.5u %3d%% ", iter, percent);
        Rprintf("Iteration: %.5u  %3d%% ", iter, percent);
        for(i = 0 ; i < nitem ; i++ ) {
          Rprintf("% .3f ", oldbeta(i));
        }
        Rprintf(" %.3f ", oldgamma);
        Rprintf(" %.3f\n", pr_sd_theta);
      }
    }else{
      // progress bar
      progressbar(iter+1,niter);
    }

  } //for end

  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["z"] = samp_z;
  output["w"] = samp_w;
  output["gamma"] = samp_gamma;
  output["pi"] = samp_pi;
  output["xi"] = samp_xi;
  output["sigma_theta"] = samp_sd_theta;
  output["sigma"] = samp_sd;
  output["map"] = samp_mle;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  output["accept_z"] = accept_z;
  output["accept_w"] = accept_w;
  output["accept_gamma"] = accept_gamma;

  return(output);

} // function end


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.

