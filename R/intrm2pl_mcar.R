#' 2pl LSIRM model using multiplicative effect for missing completely at random data. 
#' 
#' @description \link{intrm2pl_mcar} is used to fit 2pl LSIRM model using multiplicative effect in incomplete data under missing completely at random assumption. 
#' \link{intrm2pl_mcar} factorizes item response matrix into column-wise item effect, row-wise respondent effect and further embeds multiplicative effect in a latent space, while ignoring the missing element under the assumption of missing completely at random. Unlike 1pl model, 2pl model assumes the item effect can vary according to respondent, allowing additional parameter multiplied with respondent effect.  The resulting latent space provides an interaction map that represents interactions between respondents and items. 
#' 
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param ndim Numeric; dimension of latent space. default value is 2.
#' @param niter Numeric; number of iterations to run MCMC sampling. default value is 15000.
#' @param nburn Numeric; number of initial, pre-thinning, MCMC iterations to discard. default value is 2500.
#' @param nthin Numeric;number of thinning, MCMC iterations to discard. default value is 5.
#' @param nprint Numeric; MCMC samples is displayed during execution of MCMC chain for each \code{nprint}. default value is 500.
#' @param jump_beta Numeric; jumping rule of the proposal density for beta. default value is 0.4.
#' @param jump_theta Numeric; jumping rule of the proposal density for theta. default value is 1.0.
#' @param jump_alpha Numeric; jumping rule of the proposal density for alpha. default value is 1.0.
#' @param jump_delta Numeric; jumping rule of the proposal density for delta default value is 1.0.
#' @param pr_mean_beta Numeric; mean of normal prior for beta. default value is 0.
#' @param pr_sd_beta Numeric; standard deviation of normal prior for beta. default value is 1.0.
#' @param pr_mean_theta Numeric; mean of normal prior for theta. default value is 0.
#' @param pr_mean_delta Numeric; mean of normal prior for delta. default value is 0.
#' @param pr_sd_delta Numeric; standard deviation of normal prior for delta. default value is 1.0.
#' @param pr_mean_alpha Numeric; mean of normal prior for alpha. default value is 0.5.
#' @param pr_sd_alpha Numeric; mean of normal prior for beta. default value is 1.0.
#' @param pr_a_theta Numeric; shape parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param pr_b_theta Numeric; scale parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param missing Numeric; a number to replace missing values. default value is 99.
#' 
#' @return \code{intrm2pl_mcar} returns an object of  list containing the following components:
#'  \item{beta_estimate}{posterior estimation of beta.}
#'  \item{theta_estimate}{posterior estimation of theta.}
#'  \item{sigma_theta_estimate}{posterior estimation of standard deviation of theta.}
#'  \item{delta_estimate}{posterior estimation of delta.}
#'  \item{alpha_estimate}{posterior estimation of alpha.}
#'  \item{beta}{posterior samples of beta.}
#'  \item{theta}{posterior samples of theta.}
#'  \item{theta_sd}{posterior samples of standard deviation of theta.}
#'  \item{delta}{posterior samples of delta.}
#'  \item{alpha}{posterior samples of alpha.}
#'  \item{accept_beta}{accept ratio of beta.}
#'  \item{accept_theta}{accept ratio of theta.}
#'  \item{accept_alpha}{accept ratio of alpha.}
#'  \item{ls_mean_item}{posterior estimation of latent position of item.}
#'  \item{ls_mean_respondent}{posterior estimation of latent position of respondent.}
#'  \item{ls_mean_lambda}{posterior estimation lambda. The  singular value of the decomposition.}
#'  \item{ls_respondent}{posterior samples of latent positon of respondent.}
#'  \item{ls_item}{posterior samples of latent positon of item.}
#'  \item{ls_lambda}{posterior samples of lambda which is singular value of decomposition.}
#' 
#' @details \code{intrm2pl_mcar} models the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j}  in the shared metric space. For 2pl model, the the item effect is assumed to have additional discrimination parameter \eqn{\alpha_i} multiplied by \eqn{\theta_j}: \deqn{logit(P(Y_{j,i} = 1|\theta_j,\alpha_i,\beta_i,\delta_{j,i}))=\theta_j*\alpha_i+\beta_i+\delta_{j,i}} The final latent positions of respondents and items are the singular vectors of matrix with its \eqn{j,i} element \eqn{\delta_{j,i}}. Under the assumption of missing completely at random, the model ignores the missing element in doing inference. For the details of missing completely at random assumption and data augmentation, see References.
#' @references  Little, R. J., & Rubin, D. B. (2019). Statistical analysis with missing data (Vol. 793). John Wiley & Sons.
#' @examples 
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5), ncol=10, nrow=50)
#' 
#' # generate example missing indicator matrix
#' missing_mat     <- matrix(rbinom(500, size = 1, prob = 0.2), ncol=10, nrow=50)
#' 
#' # make missing value with missing indicator matrix
#' data[missing_mat==1] <- 99 
#' 
#' lsirm_result <- intrm2pl_mcar(data)
#' 
#' # The code following can achieve the same result.
#' lsirm_result <- intrm2pl(data, missing_data = 'mcar')
#' }
#' 
#' @export
intrm2pl_mcar = function(data, ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                    jump_beta = 0.4, jump_theta = 1, jump_alpha = 1.0, jump_delta = 1, 
                    pr_mean_beta = 0, pr_sd_beta = 1, pr_mean_theta = 0, 
                    pr_mean_delta = 0, pr_sd_delta = 1, pr_mean_alpha = 0.5, pr_sd_alpha = 1, 
                    pr_a_theta = 0.001, pr_b_theta = 0.001, missing = 99){
  
  output <- intrm2pl_mcar_cpp(data, niter, nburn, nthin, nprint,
                         jump_beta, jump_theta,jump_alpha,jump_delta,
                         pr_mean_beta, pr_sd_beta, pr_mean_theta, 
                         pr_mean_delta, pr_sd_delta, pr_mean_alpha, pr_sd_alpha, 
                         pr_a_theta, pr_b_theta, missing)
  
  nsample <- nrow(data)
  nitem <- ncol(data)

  beta.estimate = apply(output$beta, 2, mean)
  theta.estimate = apply(output$theta, 2, mean)
  alpha.estimate = apply(output$alpha, 2, mean)
  sigma_theta.estimate = mean(output$sigma_theta)
  delta.estimate = matrix(NA,nsample,nitem)
  
  for(k in 1:nsample){
    for(i in 1:nitem) delta.estimate[k,i] = mean(output$delta[,k,i])
  }
  
  #SVD
  ls.mean.delta = svd(delta.estimate, nu=ndim, nv=ndim)
  ls.mean.item  = ls.mean.delta$v
  ls.mean.sample = ls.mean.delta$u
  ls.mean.lambda = ls.mean.delta$d
  
  ls.item = array(NA, dim=c(nrow(output$beta),nitem,ndim))
  ls.sample = array(NA, dim=c(nrow(output$beta),nsample,ndim))
  ls.lambda = matrix(NA,nrow(output$beta),ndim)
  for(iter in 1:nrow(output$beta)){
    temp = svd(output$delta[iter,,], nu=ndim, nv=ndim)
    ls.item[iter,,] = temp$v
    ls.sample[iter,,] = temp$u
    ls.lambda[iter,] = temp$d[1:ndim]
  }
  
  return(list(beta_estimate  = beta.estimate,
              theta_estimate = theta.estimate,
              sigma_theta_estimate    = sigma_theta.estimate,
              delta_estimate = delta.estimate,
              alpha_estimate = alpha.estimate,
              beta           = output$beta,
              theta          = output$theta,
              theta_sd       = output$sigma_theta,
              delta          = output$delta,
              alpha          = output$alpha,
              accept_beta    = output$accept_beta,
              accept_theta   = output$accept_theta,
              accept_alpha   = output$accept_alpha,
              ls_mean_item   = ls.mean.item,
              ls_mean_respondent = ls.mean.sample,
              ls_mean_lambda = ls.mean.lambda,
              ls_item        = ls.item,
              ls_respondent  = ls.sample,
              ls_lambda      = ls.lambda))
}
