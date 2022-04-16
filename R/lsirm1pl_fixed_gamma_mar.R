#' 1pl LSIRM model fixing gamma to 1 for missing at random data. 
#' 
#' @description \link{lsirm1pl_fixed_gamma_mar} is used to fit LSIRM model with gamma fixed to 1 in incomplete data assumed to be missing at random. 
#' \link{lsirm1pl_fixed_gamma_mar} factorizes item response matrix into column-wise item effect, row-wise respondent effect and further embeds interaction effect in a latent space, while considering the missing element under the assumption of missing at random. The resulting latent space provides an interaction map that represents interactions between respondents and items. 
#' 
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param ndim Numeric; dimension of latent space. default value is 2.
#' @param niter Numeric; number of iterations to run MCMC sampling. default value is 15000.
#' @param nburn Numeric; number of initial, pre-thinning, MCMC iterations to discard. default value is 2500.
#' @param nthin Numeric; number of thinning, MCMC iterations to discard. default value is 5.
#' @param nprint Numeric; MCMC samples is displayed during execution of MCMC chain for each \code{nprint}. default value is 500.
#' @param jump_beta Numeric; jumping rule of the proposal density for beta. default value is 0.4.
#' @param jump_theta Numeric; jumping rule of the proposal density for theta. default value is 1.0.
#' @param jump_z Numeric; jumping rule of the proposal density for z. default value is 0.5.
#' @param jump_w Numeric; jumping rule of the proposal density for w. default value is 0.5.
#' @param pr_mean_beta Numeric; mean of normal prior for beta. default value is 0.
#' @param pr_sd_beta Numeric; standard deviation of normal prior for beta. default value is 1.0
#' @param pr_mean_theta Numeric; mean of normal prior for theta. default value is 0.
#' @param pr_a_theta Numeric; shape parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param pr_b_theta Numeric; scale parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param missing Numeric; a number to replace missing values. default value is 99.
#' 
#' @return \code{lsirm1pl_fixed_gamma_mar} returns an object of  list containing the following components:
#'  \item{beta_estimate}{posterior estimation of beta.}
#'  \item{theta_estimate}{posterior estimation of theta.}
#'  \item{sigma_theta_estimate}{posterior estimation of standard deviation of theta.}
#'  \item{z_estimate}{posterior estimation of z.}
#'  \item{w_estimate}{posterior estimation of w.}
#'  \item{imp_estimate}{probability of imputating a missing value with 1.}
#'  \item{beta}{posterior samples of beta.}
#'  \item{theta}{posterior samples of theta.}
#'  \item{theta_sd}{posterior samples of standard deviation of theta.}
#'  \item{z}{posterior samples of z. The output is 3-dimensional matrix with last axis represent the dimension of latent space.}
#'  \item{w}{posterior samples of w. The output is 3-dimensional matrix with last axis represent the dimension of latent space.}
#'  \item{imp}{imputation for missing Values using posterior samples.}
#'  \item{accept_beta}{accept ratio of beta.}
#'  \item{accept_theta}{accept ratio of theta.}
#'  \item{accept_w}{accept ratio of w.}
#'  \item{accept_z}{accept ratio of z.}
#' 
#' @details \code{lsirm1pl_fixed_gamma_mar} models the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space: \deqn{logit(P(Y_{j,i} = 1|\theta_j,\beta_i,z_j,w_i))=\theta_j+\beta_i-||z_j-w_i||} Under the assumption of missing at random, the model takes the missing element into consideration in the sampling procedure. For the details of missing at random assumption and data augmentation, see References.
#' 
#' @references Little, R. J., & Rubin, D. B. (2019). Statistical analysis with missing data (Vol. 793). John Wiley & Sons.
#' 
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
#' lsirm_result <- lsirm1pl_fixed_gamma_mar(data)
#' 
#' # The code following can achieve the same result.
#' lsirm_result <- lsirm1pl(data, spikenslab = FALSE, fixed_gamma = TRUE, missing_data = 'mar')
#' }
#' 
#' @export
lsirm1pl_fixed_gamma_mar = function(data, ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                               jump_beta = 0.4, jump_theta = 1.0, jump_z = 0.5, jump_w = 0.5,
                               pr_mean_beta = 0, pr_sd_beta = 1.0, pr_mean_theta = 0, pr_a_theta = 0.001, pr_b_theta = 0.001, missing = 99){
  
  output <- lsirm1pl_fixed_gamma_mar_cpp(data, ndim, niter, nburn, nthin, nprint,
                           jump_beta, jump_theta, jump_z, jump_w, 
                           pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_a_theta, pr_b_theta,missing)
  
  nsample <- nrow(data)
  nitem <- ncol(data)
  
  nmcmc = as.integer((niter - nburn) / nthin)
  max.address = min(which.max(output$map))
  w.star = output$w[max.address,,]
  z.star = output$z[max.address,,]
  w.proc = array(0,dim=c(nmcmc,nitem,ndim))
  z.proc = array(0,dim=c(nmcmc,nsample,ndim))

  for(iter in 1:nmcmc){
    z.iter = output$z[iter,,]
    if(iter != max.address) z.proc[iter,,] = procrustes(z.iter,z.star)$X.new
    else z.proc[iter,,] = z.iter
    
    w.iter = output$w[iter,,]
    if(iter != max.address) w.proc[iter,,] = procrustes(w.iter,w.star)$X.new
    else w.proc[iter,,] = w.iter
  }
  
  w.est = colMeans(w.proc, dims = 1)
  z.est = colMeans(z.proc, dims = 1)
  
  beta.estimate = apply(output$beta, 2, mean)
  theta.estimate = apply(output$theta, 2, mean)
  sigma_theta.estimate = mean(output$sigma_theta)
  
  return(list(beta_estimate  = beta.estimate,
              theta_estimate = theta.estimate,
              sigma_theta_estimate    = sigma_theta.estimate,
              z_estimate     = z.est,
              w_estimate     = w.est,
              beta           = output$beta,
              theta          = output$theta,
              theta_sd       = output$sigma_theta,
              z              = z.proc,
              w              = w.proc,
              accept_beta    = output$accept_beta,
              accept_theta   = output$accept_theta,
              accept_w       = output$accept_w,
              accept_z       = output$accept_z))
}