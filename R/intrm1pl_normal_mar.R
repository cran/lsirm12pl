#' 1pl LSIRM model with normal likelihood using multiplicative effect for missing at random data. 
#' 
#' @description \link{intrm1pl_normal_mar} is used to fit LSIRM model for continuous variable with 1pl using multiplicative effect in incomplete data under missing at random assumption.
#' \link{intrm1pl_normal_mar} factorizes continuous item response matrix into column-wise item effect, row-wise respondent effect and further embeds multiplicative effect in a latent space, while considering the missing element under the assumption of missing at random. The resulting latent space provides an interaction map that represents interactions between respondents and items.
#' 
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param ndim Numeric; dimension of latent space. default value is 2.
#' @param niter Numeric; number of iterations to run MCMC sampling. default value is 15000.
#' @param nburn Numeric; number of initial, pre-thinning, MCMC iterations to discard. default value is 2500.
#' @param nthin Numeric;number of thinning, MCMC iterations to discard. default value is 5.
#' @param nprint Numeric; MCMC samples is displayed during execution of MCMC chain for each \code{nprint}. default value is 500.
#' @param jump_beta Numeric; jumping rule of the proposal density for beta. default value is 0.4.
#' @param jump_theta Numeric; jumping rule of the proposal density for theta. default value is 1.0.
#' @param jump_delta Numeric; jumping rule of the proposal density for delta default value is 1.0.
#' @param pr_mean_beta Numeric; mean of normal prior for beta. default value is 0.
#' @param pr_sd_beta Numeric; standard deviation of normal prior for beta. default value is 1.0.
#' @param pr_mean_theta Numeric; mean of normal prior for theta. default value is 0.
#' @param pr_mean_delta Numeric; mean of normal prior for delta. default value is 0.
#' @param pr_sd_delta Numeric; standard deviation of normal prior for delta. default value is 1.0.
#' @param pr_a_theta Numeric; shape parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param pr_b_theta Numeric; scale parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param pr_a_eps Numeric; shape parameter of inverse gamma prior for variance of data likelihood. default value is 0.001.
#' @param pr_b_eps Numeric; scale parameter of inverse gamma prior for variance of data likelihood default value is 0.001.
#' @param missing Numeric; a number to replace missing values. default value is 99.
#' 
#' @return \code{intrm1pl_normal_mar} returns an object of  list containing the following components:
#'  \item{beta_estimate}{posterior estimation of beta.}
#'  \item{theta_estimate}{posterior estimation of theta.}
#'  \item{sigma_theta_estimate}{posterior estimation of standard deviation of theta.}
#'  \item{sigma_estimate}{posterior estimation of standard deviation.}
#'  \item{delta_estimate}{posterior estimation of delta.}
#'  \item{imp_estimate}{estimation of imputing missing values.}
#'  \item{beta}{posterior samples of beta.}
#'  \item{theta}{posterior samples of theta.}
#'  \item{theta_sd}{posterior samples of standard deviation of theta.}
#'  \item{sigma}{posterior samples of standard deviation.}
#'  \item{delta}{posterior samples of delta.}
#'  \item{imp}{imputation for missing Values using posterior samples.}
#'  \item{accept_beta}{accept ratio of beta.}
#'  \item{accept_theta}{accept ratio of theta.}
#'  \item{ls_mean_item}{posterior estimation of latent position of item.}
#'  \item{ls_mean_respondent}{posterior estimation of latent position of respondent.}
#'  \item{ls_mean_lambda}{posterior estimation lambda. The  singular value of the decomposition.}
#'  \item{ls_respondent}{posterior samples of latent positon of respondent.}
#'  \item{ls_item}{posterior samples of latent positon of item.}
#'  \item{ls_lambda}{posterior samples of lambda which is singular value of decomposition.}
#'  
#' @details \code{intrm1pl_normal_mar} models the continuous value of response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} in the shared metric space: \deqn{Y_{j,i} = \theta_j+\beta_i+\delta_{j,i} + e_{j,i}} where the error \eqn{e_{j,i} \sim N(0,\sigma^2)}. The final latent positions of respondents and items are the singular vectors of matrix with its \eqn{j,i} element \eqn{\delta_{j,i}}. Under the assumption of missing at random, the model takes the missing element into consideration in the sampling procedure. For the details of missing at random assumption and data augmentation, see References.
#' @references  Little, R. J., & Rubin, D. B. (2019). Statistical analysis with missing data (Vol. 793). John Wiley & Sons.
#'  
#' @examples 
#' # generate example (continuous) item response matrix
#' data     <- matrix(rnorm(500, mean = 0, sd = 1),ncol=10,nrow=50)
#' 
#' # generate example missing indicator matrix
#' missing_mat     <- matrix(rbinom(500, size = 1, prob = 0.2),ncol=10,nrow=50)
#' 
#' # make missing value with missing indicator matrix
#' data[missing_mat==1] <- 99 
#' 
#' lsirm_result <- intrm1pl_normal_mar(data)
#' 
#' # The code following can achieve the same result.
#' lsirm_result <- intrm1pl_normal(data, missing_data = 'mar')
#' 
#' @export
intrm1pl_normal_mar = function(data, ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                                jump_beta = 0.4, jump_theta = 1, jump_delta = 1, 
                                pr_mean_beta = 0, pr_sd_beta = 1, pr_mean_theta = 0, 
                                pr_mean_delta = 0, pr_sd_delta = 1,
                                pr_a_theta = 0.001, pr_b_theta = 0.001, pr_a_eps = 0.001, pr_b_eps = 0.001,
                                missing = 99){
  
  output <-  intrm1pl_normal_mar_cpp(data, niter, nburn, nthin, nprint,
                                      jump_beta, jump_theta, jump_delta,
                                      pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_mean_delta, pr_sd_delta,
                                      pr_a_eps, pr_b_eps, pr_a_theta, pr_b_theta, missing)
  
  beta.estimate = apply(output$beta, 2, mean)
  theta.estimate = apply(output$theta, 2, mean)
  sigma_theta.estimate = mean(output$sigma_theta)
  sigma.estimate = mean(output$sigma)
  delta.estimate = matrix(NA,nrow(data),ncol(data))
  imp.estimate = apply(output$impute, 2, mean)
  nsample <- nrow(data)
  nitem <- ncol(data)
  
  
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
              sigma_estimate    = sigma.estimate,
              delta_estimate = delta.estimate,
              imp_estimate   = imp.estimate,
              beta           = output$beta,
              theta          = output$theta,
              theta_sd       = output$sigma_theta,
              sigma          = output$sigma,
              delta          = output$delta,
              imp            = output$impute,              
              accept_beta    = output$accept_beta,
              accept_theta   = output$accept_theta,
              ls_mean_item   = ls.mean.item,
              ls_mean_respondent = ls.mean.sample,
              ls_mean_lambda = ls.mean.lambda,
              ls_item        = ls.item,
              ls_respondent  = ls.sample,
              ls_lambda      = ls.lambda))
}
