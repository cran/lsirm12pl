#' Summary the result of LSIRM model
#'
#' @description \link{summary} is used to summary the result of LSIRM model.
#'
#' @param x object of class \code{lsirm1pl}, \code{lsirm2pl}.
#' @param \dots \dots{}
#'
#' @return \code{plot_latent} returns the plot of latent space visualize an interaction map that represents interactions between respondents and items.
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm1pl(data = data)
#' summary(lsirm_result)
#' }
#' @export
print.summary.lsirm <- function(x, ...){
  cat("==========================","\n")
  cat("Summary of model","\n")
  cat("==========================","\n\n")
  cat("Call:\t", deparse(x$call), "\n")
  cat("Model:\t", x$method, "\n")
  cat("Data type:\t", x$dtype, "\n")
  cat("Missing:\t", x$missing, "\n")
  cat(sprintf("MCMC sample of size %i, after burnin of %i iteration",
              x$mcmc.opt$niter,x$mcmc.opt$nburn),"\n\n")
  cat("Covariate coefficients posterior means: ","\n\n")
  printCoefmat(x$coef)
  cat("\n---------------------------","\n\n")
  cat("Overall BIC (Smaller is better) :",x$BIC,"\n")
  cat("\nMaximum Posterior Iteration: ", "\n")
  printCoefmat(x$map.inf)
}
