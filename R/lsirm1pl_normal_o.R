#' 1PL LSIRM with normal likelihood.
#'
#' @description \link{lsirm1pl_normal_o} is used to fit LSIRM for continuous variable with 1pl.
#' \link{lsirm1pl_normal_o} factorizes continuous item response matrix into column-wise item effect, row-wise respondent effect and further embeds interaction effect in a latent space. The resulting latent space provides an interaction map that represents interactions between respondents and items.
#'
#' @inheritParams lsirm1pl
#' @param jump_gamma Numeric; the jumping rule for the gamma proposal density. Default is 0.025
#' @param pr_mean_gamma Numeric; mean of log normal prior for gamma. Default is 0.5.
#' @param pr_sd_gamma Numeric; standard deviation of log normal prior for gamma. Default is 1.0.
#' @param pr_a_eps Numeric; the shape parameter of inverse gamma prior for variance of data likelihood. Default is 0.001.
#' @param pr_b_eps Numeric; the scale parameter of inverse gamma prior for variance of data likelihood. Default is 0.001.
#' @param verbose Logical; If TRUE, MCMC samples are printed for each \code{nprint}. Default is FALSE.
#'
#'
#' @return \code{lsirm1pl_normal_o} returns an object of  list containing the following components:
#' \item{data}{Data frame or matrix containing the variables used in the model.}
#' \item{bic}{A numeric value representing the Bayesian Information Criterion (BIC).}
#' \item{mcmc_inf}{Details about the number of MCMC iterations, burn-in periods, and thinning intervals.}
#' \item{map_inf}{The log maximum a posteriori (MAP) value and the iteration number at which this MAP value occurs.}
#' \item{beta_estimate}{Posterior estimates of the beta parameter.}
#' \item{theta_estimate}{Posterior estimates of the theta parameter.}
#' \item{sigma_theta_estimate}{Posterior estimates of the standard deviation of theta.}
#' \item{gamma_estimate}{Posterior estimates of gamma parameter.}
#' \item{z_estimate}{Posterior estimates of the z parameter.}
#' \item{w_estimate}{Posterior estimates of the w parameter.}
#' \item{beta}{Posterior samples of the beta parameter.}
#' \item{theta}{Posterior samples of the theta parameter.}
#' \item{theta_sd}{Posterior samples of the standard deviation of theta.}
#' \item{gamma}{Posterior samples of the gamma parameter.}
#' \item{z}{Posterior samples of the z parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{w}{Posterior samples of the w parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{accept_beta}{Acceptance ratio for the beta parameter.}
#' \item{accept_theta}{Acceptance ratio for the theta parameter.}
#' \item{accept_z}{Acceptance ratio for the z parameter.}
#' \item{accept_w}{Acceptance ratio for the w parameter.}
#' \item{accept_gamma}{Acceptance ratio for the gamma parameter.}
#'  \item{sigma_estimate}{Posterior estimates of the standard deviation.}
#'  \item{sigma}{Posterior samples of the standard deviation.}
#'
#' @details \code{lsirm1pl_normal_o} models the continuous value of response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term: \deqn{Y_{j,i} = \theta_j+\beta_i-\gamma||z_j-w_i|| + e_{j,i}} where the error \eqn{e_{j,i} \sim N(0,\sigma^2)}.
#' @examples
#' # generate example (continuous) item response matrix
#' data     <- matrix(rnorm(500, mean = 0, sd = 1),ncol=10,nrow=50)
#'
#' lsirm_result <- lsirm1pl_normal_o(data)
#'
#' # The code following can achieve the same result.
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = FALSE, fixed_gamma = FALSE))
#'
#' @export
lsirm1pl_normal_o = function(data, ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                             jump_beta = 0.4, jump_theta = 1.0, jump_gamma = 1.0, jump_z = 0.5, jump_w = 0.5,
                             pr_mean_beta = 0, pr_sd_beta = 1.0, pr_mean_theta = 0, pr_mean_gamma = 0.5, pr_sd_gamma =1.0,
                             pr_a_theta = 0.001, pr_b_theta = 0.001, pr_a_eps = 0.001, pr_b_eps = 0.001, verbose=FALSE){
  if(niter < nburn){
    stop("niter must be greater than burn-in process.")
  }
  if(is.data.frame(data)){
    cname = colnames(data)
  }else{
    cname = paste("item", 1:ncol(data), sep=" ")
  }

  # cat("\n\nFitting with MCMC algorithm\n")


  output <- lsirm1pl_normal_cpp(as.matrix(data), ndim, niter, nburn, nthin, nprint,
                                jump_beta, jump_theta, jump_gamma, jump_z, jump_w,
                                pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_mean_gamma, pr_sd_gamma,
                                pr_a_theta, pr_b_theta, pr_a_eps, pr_b_eps, verbose=verbose)

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

  beta.summary = data.frame(cbind(apply(output$beta, 2, mean), t(apply(output$beta, 2, function(x) quantile(x, probs = c(0.025, 0.975))))))
  colnames(beta.summary) <- c("Estimate", "2.5%", "97.5%")
  rownames(beta.summary) <- cname

  # Calculate BIC
  # cat("\n\nCalculate BIC\n")
  log_like = log_likelihood_normal_cpp(as.matrix(data), ndim, as.matrix(beta.estimate), as.matrix(theta.estimate), gamma.estimate, z.est, w.est, sigma.estimate, 99)
  p = nitem + nsample + 1 + 1 + ndim * nitem + ndim * nsample + 1
  bic = -2 * log_like[[1]] + p * log(nsample * nsample)

  result <- list(data = data,
              bic = bic,
                 mcmc_inf = mcmc.inf,
                 map_inf = map.inf,
                 beta_estimate  = beta.estimate,
                 beta_summary = beta.summary,
                 theta_estimate = theta.estimate,
                 sigma_theta_estimate    = sigma_theta.estimate,
                 gamma_estimate = gamma.estimate,
                 sigma_estimate = sigma.estimate,
                 z_estimate     = z.est,
                 w_estimate     = w.est,
                 beta           = output$beta,
                 theta          = output$theta,
                 theta_sd       = output$sigma_theta,
                 sigma          = output$sigma,
                 gamma          = output$gamma,
                 z              = z.proc,
                 w              = w.proc,
                 accept_beta    = output$accept_beta,
                 accept_theta   = output$accept_theta,
                 accept_w       = output$accept_w,
                 accept_z       = output$accept_z,
                 accept_gamma   = output$accept_gamma)
  class(result) = "lsirm"

  return(result)
}
