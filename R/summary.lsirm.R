#' Summary the result of LSIRM
#'
#' @description \link{summary} is used to summary the result of LSIRM.
#'
#' @param object Object of class \code{lsirm}.
#' @param chain.idx Numeric; Index of MCMC chain. Default is 1.
#' @param estimate Character; Specifies the type of posterior estimate to provide for beta parameters. Options are \code{"mean"}, \code{"median"}, or \code{"mode"}. Default is \code{"mean"}.
#' @param CI Numeric; The significance level for the highest posterior density interval (HPD) for the beta parameters. Default is 0.95.
#' @param \dots Additional arguments.
#'
#' @return \code{summary.lsirm} contains following elements. A print method is available.
#' \item{call}{R call used to fit the model.}
#' \item{coef}{Covariate coefficients posterior means.}
#' \item{mcmc.opt}{The number of mcmc iteration, burn-in periods, and thinning intervals.}
#' \item{map.inf}{Value of log maximum a posterior and iteration number which have log maximum a posterior.}
#' \item{BIC}{Numeric value with the corresponding Bayesian information criterion (BIC).}
#' \item{method}{Which model is fitted.}
#' \item{missing}{The assumed missing type. One of NA, "mar" and "mcar".}
#' \item{dtype}{Type of input data (Binary or Continuous).}
#' \item{ss}{Whether a model selection approach using the spike-slab prior is applied.}
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#'
#' # 1PL LSIRM object
#' lsirm_result <- lsirm(data ~ lsirm1pl())
#' summary(lsirm_result)
#' }
#' @rdname summary.lsirm
#' @export
summary.lsirm <- function(object, chain.idx = 1, estimate = 'mean', CI = 0.95, ...)
{
  if(object$method == "lsirm1pl") method = "lpl LSIRM"
  if(object$method == "lsirm2pl") method = "2pl LSIRM"
  if(object$method == "lsirmgrm") method = "ordinal GRM LSIRM"
  if(object$method == "lsirmgrm2pl") method = "ordinal GRM 2PL LSIRM"

  if(object$chains == 1){

    if(object$method %in% c("lsirmgrm", "lsirmgrm2pl")){
      # For GRM modles, we use the flattened beta_matrix for summaries
      beta_samples <- if(!is.null(object$beta_matrix)) object$beta_matrix else object$beta
      cname <- colnames(beta_samples)
      if(is.null(cname)) cname <- paste0("beta", 1:ncol(beta_samples))
    }else{
      beta_samples <- object$beta
      if(is.data.frame(object$data)){
        cname = colnames(object$data)
      }else{
        cname = paste("item", 1:ncol(object$data), sep=" ")
      }
    }

    if(estimate %in% c('mean', 'median', 'mode')){
      if(estimate == 'mean'){
        est = apply(beta_samples, 2, mean)
      }else if(estimate == 'median'){
        est = apply(beta_samples, 2, median)
      }else if(estimate == 'mode'){
        # Basic mode calculation: most frequent value (might not be ideal for continuous)
        # For now, let's use mean as a fallback or if we want true mode we'd need a density approach
        # Given this is MCMC, mean or median are usually preferred.
        # But if the user asked for mode, we should try to provide it or warn.
        est = apply(beta_samples, 2, function(x) {
          d <- density(x)
          d$x[which.max(d$y)]
        })
      }
    }else{
      warning("Estimate type is not 'mean', 'median', or 'mode'. Therefore, calculating mean estimate.")
      estimate <- 'mean'
      est = apply(beta_samples, 2, mean)
    }

    if(length(CI) == 1){
      if (CI < 0.5) {
        warning("CI should be greater than 0.5; defaulting to 0.95.")
        CI <- 0.95  # 기본값으로 재설정
        quant <- t(apply(beta_samples, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
        ci.temp = c(0.025, 0.975)
      }else{
        quant <- t(apply(beta_samples, 2, function(x) quantile(x, probs = c((1 - CI)/2, 1-(1 - CI)/2))))
        ci.temp = c((1-CI)/2, (1+CI)/2)
      }
    }else if(length(CI) == 2){
      if(CI[1]<CI[2]){
        quant <- t(apply(beta_samples, 2, function(x) quantile(x, probs = c(CI[1],CI[2]))))
        ci.temp = CI
      }else{
        stop("The lower CI bound must be less than the upper CI bound.")
      }
    }
    beta.summary = data.frame(cbind(est,quant))

    colnames(beta.summary) <- c(paste0("Estimate.", estimate),
                                paste0(ci.temp[1]*100, '%'),
                                paste0(ci.temp[2]*100, '%'))

    rownames(beta.summary) <- cname

    res <- list(call = object$call,
                coef = beta.summary,
                mcmc.opt = object$mcmc_inf,
                map.inf = object$map_inf,
                BIC = object$bic,
                method = method,
                missing = object$missing,
                dtype = object$dtype,
                ss = object$varselect,
                tuning = object$tuning,
                n.chains = 1)
  }else{

    object.chain = object[[chain.idx]]
    if(object.chain$method %in% c("lsirmgrm", "lsirmgrm2pl")){
      beta_samples <- if(!is.null(object.chain$beta_matrix)) object.chain$beta_matrix else object.chain$beta
      cname <- colnames(beta_samples)
      if(is.null(cname)) cname <- paste0("beta", 1:ncol(beta_samples))
    }else{
      beta_samples <- object.chain$beta
      if(is.data.frame(object.chain$data)){
        cname = colnames(object.chain$data)
      }else{
        cname = paste("item", 1:ncol(object.chain$data), sep=" ")
      }
    }

    if(estimate %in% c('mean', 'median', 'mode')){
      if(estimate == 'mean'){
        est = apply(beta_samples, 2, mean)
      }else if(estimate == 'median'){
        est = apply(beta_samples, 2, median)
      }else if(estimate == 'mode'){
        est = apply(beta_samples, 2, function(x) {
          d <- density(x)
          d$x[which.max(d$y)]
        })
      }
    }else{
      warning("Estimate type is not 'mean', 'median', or 'mode'. Therefore, calculating mean estimate.")
      estimate <- 'mean'
      est = apply(beta_samples, 2, mean)
    }

    if(length(CI) == 1){
      if (CI < 0.5) {
        warning("CI should be greater than 0.5; defaulting to 0.95.")
        CI <- 0.95  # 기본값으로 재설정
        quant <- t(apply(beta_samples, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
        ci.temp = c(0.025, 0.975)
      }else{
        quant <- t(apply(beta_samples, 2, function(x) quantile(x, probs = c((1 - CI)/2, 1-(1 - CI)/2))))
        ci.temp = c((1-CI)/2, (1+CI)/2)
      }
    }else if(length(CI) == 2){
      if(CI[1]<CI[2]){
        quant <- t(apply(beta_samples, 2, function(x) quantile(x, probs = c(CI[1],CI[2]))))
        ci.temp = CI
      }else{
        stop("The lower CI bound must be less than the upper CI bound.")
      }
    }
    beta.summary = data.frame(cbind(est,quant))

    colnames(beta.summary) <- c(paste0("Estimate.", estimate),
                                paste0(ci.temp[1]*100, '%'),
                                paste0(ci.temp[2]*100, '%'))
    rownames(beta.summary) <- cname

    res <- list(call = object.chain$call,
                coef = beta.summary,
                mcmc.opt = object.chain$mcmc_inf,
                map.inf = object.chain$map_inf,
                BIC = object.chain$bic,
                method = method,
                missing = object.chain$missing,
                dtype = object.chain$dtype,
                ss = object.chain$varselect,
                tuning = object.chain$tuning,
                n.chains = object$chains,
                chain = chain.idx)
  }

  class(res) <- "summary.lsirm"
  res
}
