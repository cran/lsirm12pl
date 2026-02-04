#' 2PL GRM LSIRM with spike-and-slab prior.
#'
#' @description \link{lsirmgrm2pl_ss} is used to fit 2PL GRM LSIRM with spike-and-slab prior.
#'
#' @inheritParams lsirmgrm2pl
#' @param verbose Logical; If TRUE, MCMC samples are printed for each \code{nprint}. Default is FALSE.
#' @param verbose Logical; If TRUE, MCMC samples are printed for each \code{nprint}. Default is FALSE.
#'
#' @examples
#' \donttest{
#' # generate example ordinal item response matrix
#' set.seed(123)
#' nsample <- 50
#' nitem <- 10
#' data <- matrix(sample(1:5, nsample * nitem, replace = TRUE), nrow = nsample)
#'
#' # Fit 2PL GRM LSIRM with Spike-and-Slab
#' fit <- lsirmgrm2pl_ss(data, niter = 1000, nburn = 500)
#' summary(fit)
#' }
#' 
#' @export
lsirmgrm2pl_ss <- function(data, ncat = NULL, missing_data = NA, missing.val = 99,
                          chains = 1, multicore = 1, seed = NA,
                          ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                          jump_beta = 0.4, jump_theta = 1, jump_alpha = 1, jump_gamma = 0.2, jump_z = 0.5, jump_w = 0.5,
                          pr_mean_beta = 0, pr_sd_beta = 1, pr_mean_theta = 0, pr_sd_theta = 1,
                          pr_mean_alpha = 0.5, pr_sd_alpha = 1,
                          pr_a_theta = 0.001, pr_b_theta = 0.001,
                          pr_spike_mean = -3, pr_spike_sd = 1,
                          pr_slab_mean = 0.5, pr_slab_sd = 1,
                          pr_xi_a = 1, pr_xi_b = 1,
                          adapt = NULL,
                          verbose = FALSE, fix_theta_sd = FALSE, fix_alpha_1 = TRUE) {
  lsirmgrm2pl(data = data, ncat = ncat, missing_data = missing_data, missing.val = missing.val,
              chains = chains, multicore = multicore, seed = seed,
              ndim = ndim, niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
              jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
              pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
              pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
              pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
              pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
              fixed_gamma = FALSE, spikenslab = TRUE,
              pr_spike_mean = pr_spike_mean, pr_spike_sd = pr_spike_sd,
              pr_slab_mean = pr_slab_mean, pr_slab_sd = pr_slab_sd,
              pr_xi_a = pr_xi_a, pr_xi_b = pr_xi_b,
              adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd, fix_alpha_1 = fix_alpha_1)
}
