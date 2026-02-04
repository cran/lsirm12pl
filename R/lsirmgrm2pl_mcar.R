#' 2PL GRM LSIRM for missing completely at random data.
#'
#' @description \link{lsirmgrm2pl_mcar} is used to fit 2PL GRM LSIRM in incomplete data assumed to be missing completely at random.
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
#' # generate missing value (MCAR)
#' data[sample(1:500, 50)] <- NA
#'
#' # Fit 2PL GRM LSIRM with MCAR
#' fit <- lsirmgrm2pl_mcar(data, niter = 1000, nburn = 500)
#' summary(fit)
#' }
#' 
#' @export
lsirmgrm2pl_mcar <- function(data, ncat = NULL, missing.val = 99,
                            chains = 1, multicore = 1, seed = NA,
                            ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                            jump_beta = 0.4, jump_theta = 1, jump_alpha = 1, jump_gamma = 0.2, jump_z = 0.5, jump_w = 0.5,
                            pr_mean_beta = 0, pr_sd_beta = 1, pr_mean_theta = 0, pr_sd_theta = 1,
                            pr_mean_alpha = 0.5, pr_sd_alpha = 1,
                            pr_mean_gamma = 0.5, pr_sd_gamma = 1, pr_a_theta = 0.001, pr_b_theta = 0.001,
                            adapt = NULL,
                            verbose = FALSE, fix_theta_sd = FALSE, fix_alpha_1 = TRUE) {
  lsirmgrm2pl(data = data, ncat = ncat, missing_data = "mcar", missing.val = missing.val,
              chains = chains, multicore = multicore, seed = seed,
              ndim = ndim, niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
              jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
              pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
              pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
              pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
              pr_mean_gamma = pr_mean_gamma, pr_sd_gamma = pr_sd_gamma,
              pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
              fixed_gamma = FALSE, spikenslab = FALSE,
              adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd, fix_alpha_1 = fix_alpha_1)
}
