#' 1PL LSIRM with model selection approach for missing at random data.
#' @description \link{lsirm1pl_mar_ss} is used to fit 1PL LSIRM with model selection approach based on spike-and-slab priors in incomplete data assumed to be missing at random.
#' \link{lsirm1pl_mar_ss} factorizes item response matrix into column-wise item effect, row-wise respondent effect and further embeds interaction effect in a latent space, while considering the missing element under the assumption of missing at random. The resulting latent space provides an interaction map that represents interactions between respondents and items.
#'
#'
#' @inheritParams lsirm1pl
#' @param jump_gamma Numeric; the jumping rule for the theta proposal density. Default is 1.0.
#' @param pr_spike_mean Numeric; the mean of spike prior for log gamma. Default is -3.
#' @param pr_spike_sd Numeric; the standard deviation of spike prior for log gamma. Default is 1.
#' @param pr_slab_mean Numeric; the mean of spike prior for log gamma. Default is 0.5.
#' @param pr_slab_sd Numeric; the standard deviation of spike prior for log gamma. Default is is 1.
#' @param pr_xi_a Numeric; the first shape parameter of beta prior for latent variable xi. Default is 1.
#' @param pr_xi_b Numeric; the second shape parameter of beta prior for latent variable xi. Default is 1.
#' @param missing.val Numeric; a number to replace missing values. Default is 99.
#' @param verbose Logical; If TRUE, MCMC samples are printed for each \code{nprint}. Default is FALSE.
#'
#' @return \code{lsirm1pl_mar_ss} returns an object of  list containing the following components:
#'  \item{data}{Data frame or matrix containing the variables in the model.}
#'  \item{missing.val}{A number to replace missing values.}
#'  \item{bic}{Numeric value with the corresponding BIC.}
#' \item{mcmc_inf}{Details about the number of MCMC iterations, burn-in periods, and thinning intervals.}
#' \item{map_inf}{The log maximum a posteriori (MAP) value and the iteration number at which this MAP value occurs.}
#' \item{beta_estimate}{Posterior estimates of the beta parameter.}
#' \item{theta_estimate}{Posterior estimates of the theta parameter.}
#' \item{sigma_theta_estimate}{Posterior estimates of the standard deviation of theta.}
#' \item{gamma_estimate}{posterior estimates of gamma parameter.}
#' \item{z_estimate}{Posterior estimates of the z parameter.}
#' \item{w_estimate}{Posterior estimates of the w parameter.}
#' \item{beta}{Posterior samples of the beta parameter.}
#' \item{theta}{Posterior samples of the theta parameter.}
#' \item{gamma}{Posterior samples of the gamma parameter.}
#' \item{theta_sd}{Posterior samples of the standard deviation of theta.}
#' \item{z}{Posterior samples of the z parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{w}{Posterior samples of the w parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{pi}{Posterior samples of phi which is indicator of spike and slab prior. If phi is 1, log gamma follows the slab prior, otherwise follows the spike prior. }
#'  \item{imp}{Imputation for missing Values using posterior samples.}
#' \item{accept_beta}{Acceptance ratio for the beta parameter.}
#' \item{accept_theta}{Acceptance ratio for the theta parameter.}
#' \item{accept_z}{Acceptance ratio for the z parameter.}
#' \item{accept_w}{Acceptance ratio for the w parameter.}
#' \item{accept_gamma}{Acceptance ratio for the gamma parameter.}
#'  \item{pi_estimate}{Posterior estimation of phi. inclusion probability of gamma. if estimation of phi is less than 0.5, choose Rasch model with gamma = 0, otherwise latent space model with gamma > 0. }
#'  \item{imp_estimate}{Probability of imputating a missing value with 1.}
#'
#' @details \code{lsirm1pl_mar_ss} models the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term: \deqn{logit(P(Y_{j,i} = 1 |\theta_j,\beta_i,\gamma,z_j,w_i))=\theta_j+\beta_i-\gamma||z_j-w_i||} Under the assumption of missing at random, the model takes the missing element into consideration in the sampling procedure. For the details of missing at random assumption and data augmentation, see References. \code{lsirm1pl_mar_ss} model include model selection approach based on spike-and-slab priors for log gamma. For detail of spike-and-slab priors, see References.
#'
#'
#' @references Little, R. J., & Rubin, D. B. (2019). Statistical analysis with missing data (Vol. 793). John Wiley & Sons.
#' Ishwaran, H., & Rao, J. S. (2005). Spike and slab variable selection: frequentist and Bayesian strategies. The Annals of Statistics, 33(2), 730-773.
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#'
#' # generate example missing indicator matrix
#' missing_mat     <- matrix(rbinom(500, size = 1, prob = 0.2),ncol=10,nrow=50)
#'
#' # make missing value with missing indicator matrix
#' data[missing_mat==1] <- 99
#'
#' lsirm_result <- lsirm1pl_mar_ss(data)
#'
#' # The code following can achieve the same result.
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = TRUE, fixed_gamma = FALSE,
#' missing_data = 'mar', missing = 99))
#' }
#' @export
lsirm1pl_mar_ss = function(data, ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                           jump_beta = 0.4, jump_theta = 1.0, jump_gamma = 1, jump_z = 0.5, jump_w = 0.5,
                           pr_mean_beta = 0, pr_sd_beta = 1.0, pr_mean_theta = 0,
                           pr_spike_mean = -3, pr_spike_sd = 1.0, pr_slab_mean = 0.5, pr_slab_sd = 1.0,
                           pr_a_theta = 0.001, pr_b_theta = 0.001, pr_xi_a  = 1, pr_xi_b = 1,  missing.val = 99, verbose=FALSE){
  if(niter < nburn){
    stop("niter must be greater than burn-in process.")
  }
  if(is.data.frame(data)){
    cname = colnames(data)
  }else{
    cname = paste("item", 1:ncol(data), sep=" ")
  }
  # cat("\n\nFitting with MCMC algorithm\n")

  output <- lsirm1pl_mar_ss_cpp(as.matrix(data), ndim, niter, nburn, nthin, nprint,
                                jump_beta, jump_theta, jump_gamma, jump_z, jump_w,
                                pr_mean_beta, pr_sd_beta, pr_mean_theta,
                                pr_spike_mean, pr_spike_sd, pr_slab_mean, pr_slab_sd,
                                pr_a_theta, pr_b_theta, pr_beta_a = pr_xi_a, pr_beta_b = pr_xi_b, missing.val, verbose=verbose)

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
  gamma.estimate = mean(output$gamma)
  xi.estimate = mean(output$xi)
  pi.estimate = mean(output$pi)
  imp.estimate = apply(output$impute, 2, mean)

  beta.summary = data.frame(cbind(apply(output$beta, 2, mean), t(apply(output$beta, 2, function(x) quantile(x, probs = c(0.025, 0.975))))))
  colnames(beta.summary) <- c("Estimate", "2.5%", "97.5%")
  rownames(beta.summary) <- cname

  # Calculate BIC
  # cat("\n\nCalculate BIC\n")
  missing_est = ifelse(imp.estimate > 0.5, 1, 0)
  data[data == missing.val] = missing_est
  if(pi.estimate > 0.5){
    log_like = log_likelihood_cpp(as.matrix(data), ndim, as.matrix(beta.estimate), as.matrix(theta.estimate), gamma.estimate, z.est, w.est, missing.val)
  }else{
    log_like = log_likelihood_cpp(as.matrix(data), ndim, as.matrix(beta.estimate), as.matrix(theta.estimate), 0, z.est, w.est, missing.val)
  }
  p = nitem + nsample + 1 + 1 + ndim * nitem + ndim * nsample + 2
  bic = -2 * log_like[[1]] + p * log(nsample * nsample)

  result <- list(data = data,
              missing.val = missing.val,
              bic = bic,
                 mcmc_inf = mcmc.inf,
                 map_inf = map.inf,
                 beta_estimate  = beta.estimate,
                 beta_summary = beta.summary,
                 theta_estimate = theta.estimate,
                 sigma_theta_estimate    = sigma_theta.estimate,
                 gamma_estimate = gamma.estimate,
                 z_estimate     = z.est,
                 w_estimate     = w.est,
                 imp_estimate   = imp.estimate,
                 pi_estimate    = pi.estimate,
                 beta           = output$beta,
                 theta          = output$theta,
                 theta_sd       = output$sigma_theta,
                 gamma          = output$gamma,
                 z              = z.proc,
                 w              = w.proc,
                 imp            = output$impute,
                 pi             = output$pi,
                 accept_beta    = output$accept_beta,
                 accept_theta   = output$accept_theta,
                 accept_w       = output$accept_w,
                 accept_z       = output$accept_z,
                 accept_gamma   = output$accept_gamma)
  class(result) = "lsirm"

  return(result)

}
