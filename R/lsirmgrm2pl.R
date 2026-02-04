#' Fit an ordinal LSIRM with the graded response model (2PL)
#'
#' @description
#' \code{\link{lsirmgrm2pl}} fits a two-parameter logistic (2PL) extension of the ordinal latent space item response model
#' using the graded response model (GRM) with item discrimination parameters \eqn{\alpha_i}.
#' The model factorizes the item response matrix into item thresholds, respondent ability, and item discrimination,
#' while embedding interaction effects in a latent space. The resulting latent space provides an interaction map
#' that visualizes the complex relationships between respondents and items beyond traditional IRT models.
#'
#' @references
#' De Carolis, L., Kang, I., & Jeon, M. (2025). A Latent Space Graded Response Model for Likert-Scale
#' Psychological Assessments. \emph{Multivariate Behavioral Research}. \doi{10.1080/00273171.2025.2605678}
#'
#' @inheritParams lsirmgrm
#' @param jump_beta Numeric; proposal SD for GRM thresholds. Default is 0.4. During MCMC sampling, threshold proposals are constrained to maintain the ordering \eqn{\beta_{i,1} > \beta_{i,2} > \cdots > \beta_{i,K-1}} for each item.
#' @param jump_alpha Numeric; proposal SD on log-scale for \eqn{\alpha}. Default is 1.
#' @param pr_mean_alpha Numeric; log-normal prior mean for \eqn{\alpha}. Default is 0.5.
#' @param pr_sd_alpha Numeric; log-normal prior SD for \eqn{\alpha}. Default is 1.
#' @param fix_alpha_1 Logical; if TRUE, fixes \eqn{\alpha_1 = 1}. Default is TRUE.
#' @param adapt List; optional adaptive MCMC control. If not \code{NULL}, proposal standard deviations are adapted during the burn-in period to reach a target acceptance rate and are held fixed during the main MCMC sampling.
#'   When adaptation is enabled, the reported acceptance ratios in the output (\code{accept_beta}, \code{accept_theta}, \code{accept_alpha}, etc.) are computed only from iterations after burn-in, reflecting the performance of the adapted proposal distributions.
#'   Elements of the list can include:
#'   \itemize{
#'     \item \code{use_adapt}: Logical; if \code{TRUE}, adaptive MCMC is used. Default is \code{FALSE}.
#'     \item \code{adapt_interval}: Integer; the number of iterations between each update of the proposal SDs. Default is \code{100}.
#'     \item \code{adapt_rate}: Numeric; Robbins-Monro scaling constant (c) in step size formula: adapt_rate / iteration^decay_rate. Default is \code{1.0}. Valid range: any positive value. Recommended: 0.5-2.0.
#'     \item \code{decay_rate}: Numeric; Robbins-Monro decay exponent (alpha) in step size formula. Default is \code{0.5}. Valid range: (0.5, 1]. Recommended: 0.5-0.8.
#'     \item \code{target_accept}: Numeric; target acceptance rate for scalar parameters (beta, theta, gamma, alpha). Default is \code{0.44}.
#'     \item \code{target_accept_zw}: Numeric; target acceptance rate for multi-dimensional latent positions z and w. Default is \code{0.234}.
#'     \item \code{target_accept_beta/theta/alpha/gamma}: Numeric; (optional) parameter-specific target acceptance rates to override \code{target_accept}.
#'   }

#'
#' @return An object of class \code{lsirm}. For multi-chain fits, a list where each element (\code{chain1}, \code{chain2}, etc.) is a single-chain fit of class \code{lsirm}.
#'
#' If \code{missing_data = "mar"}, the returned object additionally contains \code{imp} (MCMC draws of imputed
#' responses for each missing cell) and \code{imp_estimate} (posterior mean imputation for each missing cell).
#'
#' @details
#' \code{lsirmgrm2pl} implements the 2PL extension of the Graded Response Model (GRM) in a latent space framework.
#' Let \eqn{Y_{j,i} \in \{0,\ldots,K-1\}} be the ordered categorical response of respondent \eqn{j} to item \eqn{i}.
#' The model is defined via cumulative logits:
#' \deqn{\Pr(Y_{j,i} \ge k | \theta_j, \alpha_i, \beta_{i,k}, \gamma, z_j, w_i) = \text{logit}^{-1}(\alpha_i \theta_j + \beta_{i,k} - \gamma\,\|z_j-w_i\|)}
#' for \eqn{k=1,\ldots,K-1}, where \eqn{\alpha_i} is the item discrimination parameter and \eqn{\beta_{i,k}} are item-specific thresholds that satisfy the ordering constraint \eqn{\beta_{i,1} > \beta_{i,2} > \cdots > \beta_{i,K-1}} for identifiability.
#'
#' Missing data handling:
#' \itemize{
#'   \item \code{"mcar"}: missing responses are excluded from the likelihood.
#'   \item \code{"mar"}: missing responses are imputed by data augmentation within the MCMC.
#' }
#'
#' @examples
#' # Generate example ordinal item response matrix
#' set.seed(123)
#' nsample <- 50
#' nitem <- 10
#' data <- matrix(sample(1:5, nsample * nitem, replace = TRUE), nrow = nsample)
#'
#' # Fit 2PL GRM LSIRM using direct function call
#' fit <- lsirmgrm2pl(data, niter = 1000, nburn = 500, nthin = 2)
#' summary(fit)
#'
#' # Fit with missing data (MAR)
#' data_mar <- data
#' data_mar[sample(1:length(data), 20)] <- NA
#' fit_mar <- lsirm(data_mar ~ lsirmgrm2pl(missing_data = "mar", niter = 1000, nburn = 500))
#'
#' # Fit with Spike-and-Slab prior for model selection
#' fit_ss <- lsirm(data ~ lsirmgrm2pl(spikenslab = TRUE, niter = 1000, nburn = 500))
#'
#' # Fit with adaptive MCMC for automatic tuning
#' fit_adapt <- lsirmgrm2pl(data, niter = 2000, nburn = 1000,
#'                          adapt = list(use_adapt = TRUE, adapt_interval = 50))
#' # Check adapted jump sizes and acceptance rates
#' cat("Final jump_alpha:", fit_adapt$jump_alpha, "\n")
#' cat("Acceptance rate (post-burnin):", fit_adapt$accept_alpha, "\n")
#'
#' @export
lsirmgrm2pl <- function(data, ncat = NULL, missing_data = NA, missing.val = 99,
                        chains = 1, multicore = 1, seed = NA,
                        ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                        jump_beta = 0.4, jump_theta = 1, jump_alpha = 1, jump_gamma = 0.2, jump_z = 0.5, jump_w = 0.5,
                        pr_mean_beta = 0, pr_sd_beta = 1, pr_mean_theta = 0, pr_sd_theta = 1,
                        pr_mean_alpha = 0.5, pr_sd_alpha = 1,
                        pr_mean_gamma = 0.5, pr_sd_gamma = 1, pr_a_theta = 0.001, pr_b_theta = 0.001,
                        fixed_gamma = FALSE,
                        spikenslab = FALSE,
                        pr_spike_mean = -3, pr_spike_sd = 1,
                        pr_slab_mean = 0.5, pr_slab_sd = 1,
                        pr_xi_a = 1, pr_xi_b = 1,
                        adapt = NULL,
                        verbose = FALSE, fix_theta_sd = FALSE, fix_alpha_1 = TRUE) {

  if(!is.na(seed)){
    set.seed(seed)
  }
  cat("\n Fitting ordinal LSIRM (GRM 2PL) with MCMC algorithm\n")

  if(isTRUE(fixed_gamma) && isTRUE(spikenslab)){
    stop("fixed_gamma and spikenslab cannot both be TRUE.")
  }

  fit_one <- function(){
    lsirmgrm2pl_o(data = data, ncat = ncat, missing_data = missing_data, missing.val = missing.val,
                 ndim = ndim, niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                 jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
                 pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                 pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                 pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                 pr_mean_gamma = pr_mean_gamma, pr_sd_gamma = pr_sd_gamma,
                 pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                 fixed_gamma = fixed_gamma,
                 spikenslab = spikenslab,
                 pr_spike_mean = pr_spike_mean, pr_spike_sd = pr_spike_sd,
                 pr_slab_mean = pr_slab_mean, pr_slab_sd = pr_slab_sd,
                 pr_xi_a = pr_xi_a, pr_xi_b = pr_xi_b,
                 adapt = adapt,
                 verbose = verbose, fix_theta_sd = fix_theta_sd, fix_alpha_1 = fix_alpha_1)
  }

  if(chains > 1){
    if(multicore <= 1){
      output <- vector("list", chains)
      for(i in 1:chains){
        output[[i]] <- fit_one()
        output[[i]]$dtype <- "ordinal"
        cat(sprintf("Chain %d / %d completed\n", i, chains))
      }
    }else{
      if(chains < multicore){
        stop("Error: The number of chains must not be less than the number of cores.")
      }
      cl <- parallel::makeCluster(multicore)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      if(!is.na(seed)) parallel::clusterSetRNGStream(cl, seed)

      output <- parallel::parLapply(cl, seq_len(chains), function(i){
        lsirmgrm2pl_o(data = data, ncat = ncat, missing_data = missing_data, missing.val = missing.val,
                     ndim = ndim, niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                     jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
                     pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                     pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                     pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                     pr_mean_gamma = pr_mean_gamma, pr_sd_gamma = pr_sd_gamma,
                     pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                     fixed_gamma = fixed_gamma,
                     spikenslab = spikenslab,
                     pr_spike_mean = pr_spike_mean, pr_spike_sd = pr_spike_sd,
                     pr_slab_mean = pr_slab_mean, pr_slab_sd = pr_slab_sd,
                     pr_xi_a = pr_xi_a, pr_xi_b = pr_xi_b,
                     adapt = adapt,
                     verbose = verbose, fix_theta_sd = fix_theta_sd, fix_alpha_1 = fix_alpha_1)
      })
      for(i in 1:chains) output[[i]]$dtype <- "ordinal"
    }

    for(i in 1:chains){
      output[[i]]$call <- match.call()
      output[[i]]$method <- "lsirmgrm2pl"
      output[[i]]$missing <- output[[i]]$missing
      output[[i]]$varselect <- isTRUE(spikenslab)
      output[[i]]$fixed_gamma <- isTRUE(fixed_gamma)
      output[[i]]$chains <- chains
      class(output[[i]]) <- "lsirm"
    }

    names(output) <- paste0("chain", 1:chains)
    output$method <- "lsirmgrm2pl"
    output$chains <- chains
    class(output) <- "lsirm"
    return(output)
  }

  output <- fit_one()
  output$dtype <- "ordinal"
  output$call <- match.call()
  output$method <- "lsirmgrm2pl"
  output$varselect <- isTRUE(spikenslab)
  output$fixed_gamma <- isTRUE(fixed_gamma)
  output$chains <- chains
  class(output) <- "lsirm"
  output
}


#' Fit a single-chain ordinal LSIRM (GRM 2PL)
#'
#' @inheritParams lsirmgrm2pl
#' @return A list containing MCMC draws and posterior summaries, including:
#' \itemize{
#'   \item \code{beta}, \code{theta}, \code{gamma}, \code{alpha}, \code{z}, \code{w}: MCMC draws.
#'   \item \code{beta_estimate}, \code{theta_estimate}, \code{gamma_estimate}, \code{alpha_estimate}, \code{z_estimate},
#'   \code{w_estimate}: posterior means.
#' }
#'
#' If \code{missing_data = "mar"}, the list additionally includes \code{imp} and \code{imp_estimate} for the imputed
#' ordinal responses.
#' @export
lsirmgrm2pl_o <- function(data, ncat = NULL, missing_data = NA, missing.val = 99,
                          ndim = 2, niter = 15000, nburn = 2500, nthin = 5, nprint = 500,
                          jump_beta = 0.4, jump_theta = 1, jump_alpha = 1, jump_gamma = 0.2, jump_z = 0.5, jump_w = 0.5,
                          pr_mean_beta = 0, pr_sd_beta = 1,
                          pr_mean_theta = 0, pr_sd_theta = 1,
                          pr_mean_alpha = 0.5, pr_sd_alpha = 1,
                          pr_mean_gamma = 0.5, pr_sd_gamma = 1,
                          pr_a_theta = 0.001, pr_b_theta = 0.001,
                          fixed_gamma = FALSE,
                          spikenslab = FALSE,
                          pr_spike_mean = -3, pr_spike_sd = 1,
                          pr_slab_mean = 0.5, pr_slab_sd = 1,
                          pr_xi_a = 1, pr_xi_b = 1,
                          adapt = NULL,
                          verbose = FALSE, fix_theta_sd = FALSE, fix_alpha_1 = TRUE){

  if(niter < nburn){
    stop("niter must be greater than burn-in process.")
  }

  if(is.data.frame(data)){
    cname <- colnames(data)
  }else{
    cname <- paste("item", 1:ncol(data), sep = " ")
  }

  x <- as.matrix(data)
  if(anyNA(x)){
    x[is.na(x)] <- missing.val
    if(is.na(missing_data)){
      missing_data <- "mcar"
    }
  }

  is_mar <- identical(missing_data, "mar")

  obs <- as.numeric(x[x != missing.val])
  if(length(obs) == 0){
    stop("No observed responses found after removing missing values.")
  }

  y_base <- if(min(obs) >= 1) 1 else 0

  if(is.null(ncat)){
    if(min(obs) >= 1) ncat <- max(obs)
    else ncat <- max(obs) + 1
  }

  if(isTRUE(fixed_gamma) && isTRUE(spikenslab)){
    stop("fixed_gamma and spikenslab cannot both be TRUE.")
  }

  if(isTRUE(spikenslab)){
    output <- if(is_mar){
      lsirmgrm2pl_ss_mar_cpp(data = x, ndim = ndim, ncat = ncat,
                             niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                             jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
                             pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                             pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                             pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                             pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                             pr_spike_mean = pr_spike_mean, pr_spike_sd = pr_spike_sd,
                             pr_slab_mean = pr_slab_mean, pr_slab_sd = pr_slab_sd,
                             pr_beta_a = pr_xi_a, pr_beta_b = pr_xi_b,
                             missing = missing.val, adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd,
                             fix_alpha_1 = fix_alpha_1)
    }else{
      lsirmgrm2pl_ss_cpp(data = x, ndim = ndim, ncat = ncat,
                         niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                         jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
                         pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                         pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                         pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                         pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                         pr_spike_mean = pr_spike_mean, pr_spike_sd = pr_spike_sd,
                         pr_slab_mean = pr_slab_mean, pr_slab_sd = pr_slab_sd,
                         pr_beta_a = pr_xi_a, pr_beta_b = pr_xi_b,
                         missing = missing.val, adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd,
                         fix_alpha_1 = fix_alpha_1)
    }
  } else if(isTRUE(fixed_gamma)){
    output <- if(is_mar){
      lsirmgrm2pl_fixed_gamma_mar_cpp(data = x, ndim = ndim, ncat = ncat,
                                      niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                                      jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_z = jump_z, jump_w = jump_w,
                                      pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                                      pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                                      pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                                      pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                                      missing = missing.val, adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd,
                                      fix_alpha_1 = fix_alpha_1)
    }else{
      lsirmgrm2pl_fixed_gamma_cpp(data = x, ndim = ndim, ncat = ncat,
                                  niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                                  jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_z = jump_z, jump_w = jump_w,
                                  pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                                  pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                                  pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                                  pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                                  missing = missing.val, adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd,
                                  fix_alpha_1 = fix_alpha_1)
    }
  } else {
    output <- if(is_mar){
      lsirmgrm2pl_mar_cpp(data = x, ndim = ndim, ncat = ncat,
                          niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                          jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
                          pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                          pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                          pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                          pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                          pr_mean_gamma = pr_mean_gamma, pr_sd_gamma = pr_sd_gamma,
                          missing = missing.val, adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd,
                          fix_alpha_1 = fix_alpha_1)
    }else{
      lsirmgrm2pl_cpp(data = x, ndim = ndim, ncat = ncat,
                      niter = niter, nburn = nburn, nthin = nthin, nprint = nprint,
                      jump_beta = jump_beta, jump_theta = jump_theta, jump_alpha = jump_alpha, jump_gamma = jump_gamma, jump_z = jump_z, jump_w = jump_w,
                      pr_mean_beta = pr_mean_beta, pr_sd_beta = pr_sd_beta,
                      pr_a_theta = pr_a_theta, pr_b_theta = pr_b_theta,
                      pr_mean_theta = pr_mean_theta, pr_sd_theta = pr_sd_theta,
                      pr_mean_alpha = pr_mean_alpha, pr_sd_alpha = pr_sd_alpha,
                      pr_mean_gamma = pr_mean_gamma, pr_sd_gamma = pr_sd_gamma,
                      missing = missing.val, adapt = adapt, verbose = verbose, fix_theta_sd = fix_theta_sd,
                      fix_alpha_1 = fix_alpha_1)
    }
  }

  mcmc.inf <- list(nburn = nburn, niter = niter, nthin = nthin)
  nsample <- nrow(x)
  nitem <- ncol(x)
  Kminus1 <- ncat - 1

  nmcmc <- as.integer((niter - nburn) / nthin)
  max.address <- min(which.max(output$map))
  map.inf <- data.frame(value = output$map[which.max(output$map)], iter = which.max(output$map))

  w.star <- output$w[max.address,,]
  z.star <- output$z[max.address,,]
  w.proc <- array(0, dim = c(nmcmc, nitem, ndim))
  z.proc <- array(0, dim = c(nmcmc, nsample, ndim))

  cat("\n")

  for(iter in 1:nmcmc){
    z.iter <- output$z[iter,,]
    w.iter <- output$w[iter,,]

    if(ndim == 1){
      z.iter <- as.matrix(z.iter)
      w.iter <- as.matrix(w.iter)
      z.star <- as.matrix(z.star)
      w.star <- as.matrix(w.star)
    }

    if(iter != max.address) z.proc[iter,,] <- MCMCpack::procrustes(z.iter, z.star)$X.new
    else z.proc[iter,,] <- z.iter

    if(iter != max.address) w.proc[iter,,] <- MCMCpack::procrustes(w.iter, w.star)$X.new
    else w.proc[iter,,] <- w.iter
  }

  w.est <- colMeans(w.proc, dims = 1)
  z.est <- colMeans(z.proc, dims = 1)

  beta.array <- output$beta
  beta.mat <- matrix(NA_real_, nrow = nmcmc, ncol = nitem * Kminus1)
  
  # Generate column names in the same order as data is stored
  beta_colnames <- character(nitem * Kminus1)
  idx <- 1
  for(i in 1:nitem){
    for(t in 1:Kminus1){
      beta_colnames[idx] <- paste0(cname[i], ":th", t)
      beta.mat[, idx] <- beta.array[, i, t]
      idx <- idx + 1
    }
  }
  colnames(beta.mat) <- beta_colnames

  beta.estimate <- apply(beta.array, c(2,3), mean)
  rownames(beta.estimate) <- cname
  colnames(beta.estimate) <- paste0("th", 1:Kminus1)

  theta.estimate <- apply(output$theta, 2, mean)
  alpha.estimate <- apply(output$alpha, 2, mean)
  sigma_theta.estimate <- mean(output$sigma_theta)
  if(isTRUE(fixed_gamma)){
    if(is.null(output$gamma)) output$gamma <- rep(1, nmcmc)
    gamma.estimate <- 1
  }else{
    gamma.estimate <- mean(output$gamma)
  }

  beta.summary <- data.frame(cbind(apply(beta.mat, 2, mean), t(apply(beta.mat, 2, function(x) quantile(x, probs = c(0.025, 0.975))))))
  colnames(beta.summary) <- c("Estimate.mean", "2.5%", "97.5%")
  rownames(beta.summary) <- colnames(beta.mat)

  miss_idx <- NULL
  imp <- NULL
  imp_estimate <- NULL
  x_bic <- x
  if(is_mar){
    miss_idx <- which(x == missing.val, arr.ind = TRUE)
    imp <- output$impute
    if(!is.null(imp) && nrow(miss_idx) > 0){
      imp_estimate <- colMeans(imp)
      imp_cat <- as.integer(round(imp_estimate))
      min_cat <- y_base
      max_cat <- y_base + ncat - 1
      imp_cat <- pmin(pmax(imp_cat, min_cat), max_cat)
      x_bic[miss_idx] <- imp_cat
    }
  }

  log_like <- log_likelihood_grm2pl_cpp(x_bic, ndim, ncat,
                                       beta_est = beta.estimate,
                                       alpha_est = as.matrix(alpha.estimate),
                                       theta_est = as.matrix(theta.estimate),
                                       gamma_est = gamma.estimate,
                                       z_est = z.est,
                                       w_est = w.est,
                                       missing = missing.val)

  p_gamma <- if(isTRUE(fixed_gamma)) 0 else 1
  p <- nitem * Kminus1 + nitem + nsample + 1 + p_gamma + ndim * nitem + ndim * nsample
  bic <- -2 * log_like[[1]] + p * log(nsample * nitem)

  if(isTRUE(fixed_gamma)){
    prior_gamma <- list(fixed = 1)
  } else if(isTRUE(spikenslab)){
    prior_gamma <- list(
      spike = c(pr_spike_mean, pr_spike_sd),
      slab = c(pr_slab_mean, pr_slab_sd),
      xi_beta = c(pr_xi_a, pr_xi_b)
    )
  } else {
    prior_gamma <- list(log_mean = pr_mean_gamma, log_sd = pr_sd_gamma)
  }

  # Convert beta matrix to list of matrices (one per threshold) for consistency with tests
  beta_list <- list()
  for(t in 1:Kminus1) {
    # Extract columns for this threshold across all items
    # The columns in beta.mat are ordered: item1:th1, item1:th2, item2:th1, item2:th2, ...
    # So for threshold t, we need columns: t, t + Kminus1, t + 2 * Kminus1, ...
    cols <- seq(t, ncol(beta.mat), by = Kminus1)
    beta_list[[t]] <- beta.mat[, cols]
    colnames(beta_list[[t]]) <- cname
  }
  names(beta_list) <- paste0("th", 1:Kminus1)

  result <- list(data = x,
                 ncat = ncat,
                 y_base = y_base,
                 missing.val = missing.val,
                 bic = bic,
                 prior = list(
                   beta = list(mean = pr_mean_beta, sd = pr_sd_beta),
                   theta = list(mean = pr_mean_theta, sd_init = pr_sd_theta, var_inv_gamma = c(pr_a_theta, pr_b_theta)),
                   alpha = list(log_mean = pr_mean_alpha, log_sd = pr_sd_alpha, fix_alpha_1 = fix_alpha_1),
                   gamma = prior_gamma,
                   z = list(mean = 0, sd = 1),
                   w = list(mean = 0, sd = 1)
                 ),
                 mcmc_inf = mcmc.inf,
                 map_inf = map.inf,
                 beta_estimate = beta.estimate,
                 beta_summary = beta.summary,
                 theta_estimate = theta.estimate,
                 alpha_estimate = alpha.estimate,
                 sigma_theta_estimate = sigma_theta.estimate,
                 gamma_estimate = gamma.estimate,
                 z_estimate = z.est,
                 w_estimate = w.est,
                 beta = beta_list,
                 beta_matrix = beta.mat,
                 beta_array = beta.array,
                 theta = output$theta,
                 theta_sd = output$sigma_theta,
                 gamma = output$gamma,
                 alpha = output$alpha,
                 z = z.proc,
                 w = w.proc,
                 z_raw = output$z,
                 w_raw = output$w,
                 accept_beta = output$accept_beta,
                 accept_theta = output$accept_theta,
                 accept_z = output$accept_z,
                 accept_w = output$accept_w,
                 accept_gamma = if(!is.null(output$accept_gamma)) output$accept_gamma else rep(NA_real_, nmcmc),
                 accept_alpha = if(!is.null(output$accept_alpha)) output$accept_alpha else rep(NA_real_, nmcmc),
                 tuning = output$tuning,
                 missing = missing_data,
                 varselect = isTRUE(spikenslab),
                 fixed_gamma = isTRUE(fixed_gamma))

  if(isTRUE(spikenslab)){
    result$pi <- as.numeric(output$pi)
    result$xi <- as.numeric(output$xi)
    result$pi_estimate <- mean(result$pi)
    result$xi_estimate <- mean(result$xi)
  }

  if(is_mar){
    result$imp_estimate <- imp_estimate
    result$imp <- imp
  }

  result
}
