#' 1PL Rasch model. 
#' 
#' @description \link{onepl} is used to fit 1PL Rasch model. 
#' 
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param niter Numeric; number of iterations to run MCMC sampling. default value is 15000.
#' @param nburn Numeric; number of initial, pre-thinning, MCMC iterations to discard. default value is 2500.
#' @param nthin Numeric;number of thinning, MCMC iterations to discard. default value is 5.
#' @param nprint Numeric; MCMC samples is displayed during execution of MCMC chain for each \code{nprint}. default value is 500.
#' @param jump_beta Numeric; jumping rule of the proposal density for beta. default value is 0.4.
#' @param jump_theta Numeric; jumping rule of the proposal density for theta. default value is 1.0.
#' @param pr_mean_beta Numeric; mean of normal prior for beta. default value is 0.
#' @param pr_sd_beta Numeric; standard deviation of normal prior for beta. default value is 1.0.
#' @param pr_mean_theta Numeric; mean of normal prior for theta. default value is 0.
#' @param pr_a_theta Numeric; shape parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' @param pr_b_theta Numeric; scale parameter of inverse gamma prior for variance of theta. default value is 0.001.
#' 
#' 
#' @return \code{onepl} returns an object of  list containing the following components:
#'  \item{beta_estimate}{posterior estimation of beta.}
#'  \item{theta_estimate}{posterior estimation of theta.}
#'  \item{sigma_theta_estimate}{posterior estimation of standard deviation of theta.}
#'  \item{beta}{posterior samples of beta.}
#'  \item{theta}{posterior samples of theta.}
#'  \item{theta_sd}{posterior samples of standard deviation of theta.}
#'  \item{accept_beta}{accept ratio of beta.}
#'  \item{accept_theta}{accept ratio of theta.}
#' 
#' @details \code{onepl} models the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j}: \deqn{logit(P(Y_{j,i} = 1|\theta_j,\beta_i))=\theta_j+\beta_i}
#' 
#' @examples 
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' 
#' result <- onepl(data)
#' }
#' @export
onepl = function(data, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                 jump_beta = 0.4, jump_theta = 1.0, 
                 pr_mean_beta = 0, pr_sd_beta = 1.0, pr_mean_theta = 0, 
                 pr_a_theta = 0.001, pr_b_theta = 0.001){
  
  if(is.data.frame(data)){
    cname = colnames(data)
  }else{
    cname = paste("item", 1:ncol(data), sep=" ")
  }
  
  output <- onepl_cpp(as.matrix(data), niter, nburn, nthin, nprint,
                      jump_beta, jump_theta, 
                      pr_mean_beta, pr_sd_beta, pr_mean_theta, 
                      pr_a_theta, pr_b_theta)
  
  mcmc.inf = list(nburn=nburn, niter=niter, nthin=nthin)
  nsample <- nrow(data)
  nitem <- ncol(data)
  
  beta.estimate = apply(output$beta, 2, mean)
  theta.estimate = apply(output$theta, 2, mean)
  sigma_theta.estimate = mean(output$sigma_theta)
  
  beta.summary = data.frame(cbind(apply(output$beta, 2, mean), t(apply(output$beta, 2, function(x) quantile(x, probs = c(0.025, 0.975))))))
  colnames(beta.summary) <- c("Estimate", "2.5%", "97.5%")
  rownames(beta.summary) <- cname
  
  result <- list(mcmc_inf = mcmc.inf,
              beta_estimate  = beta.estimate,
              beta_summary = beta.summary,
              theta_estimate = theta.estimate,
              sigma_theta_estimate    = sigma_theta.estimate,
              beta           = output$beta,
              theta          = output$theta,
              theta_sd       = output$sigma_theta,
              accept_beta    = output$accept_beta,
              accept_theta   = output$accept_theta)
  class(result) = "lsirm"
  
  return(result)            
}
