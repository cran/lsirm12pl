#' 1PL LSIRM with normal likelihood and model selection approach.
#'
#' @description \link{lsirm1pl_normal_ss} is used to fit LSIRM with model selection approach based on spike-and-slab priors for continuous variable with 1pl.
#' LSIRM factorizes continuous item response matrix into column-wise item effect, row-wise respondent effect and further embeds interaction effect in a latent space. The resulting latent space provides an interaction map that represents interactions between respondents and items.
#'
#' @inheritParams lsirm1pl
#' @param jump_gamma Numeric; the jumping rule for the theta proposal density. Default is 1.0.
#' @param pr_spike_mean Numeric; mean of spike prior for log gamma default value is -3.
#' @param pr_spike_sd Numeric; standard deviation of spike prior for log gamma default value is 1.
#' @param pr_slab_mean Numeric; mean of spike prior for log gamma default value is 0.5.
#' @param pr_slab_sd Numeric; standard deviation of spike prior for log gamma default value is 1.
#' @param pr_a_eps Numeric; shape parameter of inverse gamma prior for variance of data likelihood. Default is 0.001.
#' @param pr_b_eps Numeric; scale parameter of inverse gamma prior for variance of data likelihood. Default is 0.001.
#' @param pr_xi_a Numeric; first shape parameter of beta prior for latent variable xi. Default is 1.
#' @param pr_xi_b Numeric; second shape parameter of beta prior for latent variable xi. Default is 1.
#' @param verbose Logical; If TRUE, MCMC samples are printed for each \code{nprint}. Default is FALSE.
#'
#' @return \code{lsirm1pl_normal_ss} returns an object of  list containing the following components:
#'  \item{data}{Data frame or matrix containing the variables in the model.}
#'  \item{bic}{Numeric value with the corresponding BIC.}
#'  \item{mcmc_inf}{number of mcmc iteration, burn-in periods, and thinning intervals.}
#'  \item{map_inf}{value of log maximum a posterior and iteration number which have log maximum a posterior.}
#'  \item{beta_estimate}{posterior estimation of beta.}
#'  \item{theta_estimate}{posterior estimation of theta.}
#'  \item{sigma_theta_estimate}{posterior estimation of standard deviation of theta.}
#'  \item{sigma_estimate}{posterior estimation of standard deviation.}
#'  \item{gamma_estimate}{posterior estimation of gamma.}
#'  \item{z_estimate}{posterior estimation of z.}
#'  \item{w_estimate}{posterior estimation of w.}
#'  \item{pi_estimate}{posterior estimation of phi. inclusion probability of gamma. if estimation of phi is less than 0.5, choose Rasch model with gamma = 0, otherwise latent space model with gamma > 0. }
#'  \item{beta}{posterior samples of beta.}
#'  \item{theta}{posterior samples of theta.}
#'  \item{theta_sd}{posterior samples of standard deviation of theta.}
#'  \item{sigma}{posterior samples of standard deviation.}
#'  \item{gamma}{posterior samples of gamma.}
#'  \item{z}{posterior samples of z. The output is 3-dimensional matrix with last axis represent the dimension of latent space.}
#'  \item{w}{posterior samples of w. The output is 3-dimensional matrix with last axis represent the dimension of latent space.}
#'  \item{pi}{posterior samples of phi which is indicator of spike and slab prior. If phi is 1, log gamma follows the slab prior, otherwise follows the spike prior. }
#'  \item{accept_beta}{accept ratio of beta.}
#'  \item{accept_theta}{accept ratio of theta.}
#'  \item{accept_w}{accept ratio of w.}
#'  \item{accept_z}{accept ratio of z.}
#'  \item{accept_gamma}{accept ratio of gamma.}
#'
#' @details \code{lsirm1pl_normal_ss} models the continuous value of response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term: \deqn{Y_{j,i} = \theta_j+\beta_i-\gamma||z_j-w_i|| + e_{j,i}} where the error \eqn{e_{j,i} \sim N(0,\sigma^2)}. \code{lsrm1pl_noraml_ss} model include model selection approach based on spike-and-slab priors for log gamma. For detail of spike-and-slab priors, see References.
#' @references
#' Ishwaran, H., & Rao, J. S. (2005). Spike and slab variable selection: Frequentist and Bayesian strategies (Vol. 33). The Annals of Statistics
#'
#' @examples
#' # generate example (continuous) item response matrix
#' data     <- matrix(rnorm(500, mean = 0, sd = 1),ncol=10,nrow=50)
#'
#' lsirm_result <- lsirm1pl_normal_ss(data)
#'
#' # The code following can achieve the same result.
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = TRUE, fixed_gamma = FALSE))
#'
#' @export
lsirm1pl_normal_ss = function(data, ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                              jump_beta = 0.4, jump_theta = 1.0, jump_gamma = 1.0, jump_z = 0.5, jump_w = 0.5,
                              pr_mean_beta = 0, pr_sd_beta = 1.0, pr_mean_theta = 0,
                              pr_spike_mean = -3, pr_spike_sd = 1.0, pr_slab_mean = 0.5, pr_slab_sd = 1.0,
                              pr_a_theta = 0.001, pr_b_theta = 0.001,
                              pr_a_eps = 0.001, pr_b_eps = 0.001,
                              pr_xi_a = 0.001, pr_xi_b = 0.001, verbose=FALSE){
  if(niter < nburn){
    stop("niter must be greater than burn-in process.")
  }
  if(is.data.frame(data)){
    cname = colnames(data)
  }else{
    cname = paste("item", 1:ncol(data), sep=" ")
  }

  # cat("\n\nFitting with MCMC algorithm\n")


  output <- lsirm1pl_normal_ss_cpp(as.matrix(data), ndim, niter, nburn, nthin, nprint,
                                   jump_beta, jump_theta, jump_gamma, jump_z, jump_w,
                                   pr_mean_beta, pr_sd_beta, pr_mean_theta,
                                   pr_spike_mean, pr_spike_sd, pr_slab_mean, pr_slab_sd,
                                   pr_a_theta, pr_b_theta,
                                   pr_a_eps, pr_b_eps,
                                   pr_xi_a, pr_xi_b, verbose=verbose)

  mcmc.inf = list(nburn=nburn, niter=niter, nthin=nthin)
  nsample <- nrow(data)
  nitem <- ncol(data)

  nmcmc = as.integer((niter - nburn) / nthin)
  max.address = min(which.max(output$map))
  map.inf = data.frame(value = output$map[which.max(output$map)], iter = which.max(output$map))
  w.star = output$w[max.address,,]
  z.star = output$z[max.address,,]
  w.proc = array(0,dim=c(nmcmc,nitem,ndim))
  z.proc = array(0,dim=c(nmcmc,nsample,ndim))

  # cat("\n\nProcrustes Matching Analysis\n")
cat("\n")

  for(iter in 1:nmcmc){
    z.iter = output$z[iter,,]
    w.iter = output$w[iter,,]

    if(ndim == 1){
      z.iter = as.matrix(z.iter)
      w.iter = as.matrix(w.iter)
      z.star = as.matrix(z.star)
      w.star = as.matrix(w.star)
    }

    if(iter != max.address) z.proc[iter,,] = procrustes(z.iter,z.star)$X.new
    else z.proc[iter,,] = z.iter

    if(iter != max.address) w.proc[iter,,] = procrustes(w.iter,w.star)$X.new
    else w.proc[iter,,] = w.iter
  }

  w.est = colMeans(w.proc, dims = 1)
  z.est = colMeans(z.proc, dims = 1)

  beta.estimate = apply(output$beta, 2, mean)
  theta.estimate = apply(output$theta, 2, mean)
  sigma_theta.estimate = mean(output$sigma_theta)
  sigma.estimate = mean(output$sigma)
  gamma.estimate = mean(output$gamma)
  pi.estimate = mean(output$pi)
  xi.estimate = mean(output$xi)

  beta.summary = data.frame(cbind(apply(output$beta, 2, mean), t(apply(output$beta, 2, function(x) quantile(x, probs = c(0.025, 0.975))))))
  colnames(beta.summary) <- c("Estimate", "2.5%", "97.5%")
  rownames(beta.summary) <- cname

  # Calculate BIC
  # cat("\n\nCalculate BIC\n")
  if(pi.estimate > 0.5){
    log_like = log_likelihood_normal_cpp(as.matrix(data), ndim, as.matrix(beta.estimate), as.matrix(theta.estimate), gamma.estimate, z.est, w.est, sigma.estimate, 99)
  }else{
    log_like = log_likelihood_normal_cpp(as.matrix(data), ndim, as.matrix(beta.estimate), as.matrix(theta.estimate), 0, z.est, w.est, sigma.estimate, 99)
  }
  p = nitem + nsample + 1 + 1 + ndim * nitem + ndim * nsample + 2 + 1
  bic = -2 * log_like[[1]] + p * log(nsample * nsample)

  result <- list(data = data,
              bic = bic,
                 mcmc_inf = mcmc.inf,
                 map_inf = map.inf,
                 beta_estimate  = beta.estimate,
                 beta_summary = beta.summary,
                 theta_estimate = theta.estimate,
                 sigma_theta_estimate    = sigma_theta.estimate,
                 sigma_estimate    = sigma.estimate,
                 gamma_estimate = gamma.estimate,
                 z_estimate     = z.est,
                 w_estimate     = w.est,
                 pi_estimate    = pi.estimate,
                 beta           = output$beta,
                 theta          = output$theta,
                 theta_sd       = output$sigma_theta,
                 sigma       = output$sigma,
                 gamma          = output$gamma,
                 z              = z.proc,
                 w              = w.proc,
                 pi             = output$pi,
                 accept_beta    = output$accept_beta,
                 accept_theta   = output$accept_theta,
                 accept_w       = output$accept_w,
                 accept_z       = output$accept_z,
                 accept_gamma   = output$accept_gamma)
  class(result) = "lsirm"

  return(result)
}
