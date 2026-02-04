#' 1PL GRM LSIRM fixing gamma to 1.
#'
#' @description \link{lsirmgrm_fixed_gamma} is used to fit 1PL GRM LSIRM with gamma fixed to 1.
#' 
#' @inheritParams lsirmgrm
#' @param verbose Logical; If TRUE, MCMC samples are printed for each \code{nprint}. default value is FALSE
#'
#' @examples
#' \donttest{
#' # generate example ordinal item response matrix
#' set.seed(123)
#' nsample <- 50
#' nitem <- 10
#' data <- matrix(sample(1:5, nsample * nitem, replace = TRUE), nrow = nsample)
#'
#' # Fit 1PL GRM LSIRM with fixed gamma = 1
#' fit <- lsirmgrm_fixed_gamma(data, niter = 1000, nburn = 500)
#' summary(fit)
#' }
#' 
#' @return \code{lsirmgrm_fixed_gamma} returns an object of  list containing the same components as \code{\link{lsirmgrm}}.
#'
#' @export
lsirmgrm_fixed_gamma <- function(data, ncat = NULL, missing_data = NA, missing.val = 99,
                                chains = 1, multicore = 1, seed = NA,
                                ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                                jump_beta = 0.4, jump_theta = 1, jump_z = 0.5, jump_w = 0.5,
                                pr_mean_beta = 0, pr_sd_beta = 1, pr_mean_theta = 0, pr_sd_theta = 1,
                                pr_a_theta = 0.001, pr_b_theta = 0.001,
                                adapt = NULL,
                                verbose = FALSE, fix_theta_sd = FALSE) {
  lsirmgrm(data = data, ncat = ncat, missing_data = missing_data, missing.val = missing.val,
           chains = chains, multicore = multicore, seed = seed,
           ndim = ndim, niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
           jump_beta = jump_beta, jump_theta = jump_theta, jump_z = jump_z, jump_w = jump_w,
           pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta, 
           pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
           pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
           fixed_gamma = TRUE, spikenslab = FALSE,
           adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd)
}
