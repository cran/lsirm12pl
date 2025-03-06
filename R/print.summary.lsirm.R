#' Print the summary the result of LSIRM
#'
#' @description \link{print.summary.lsirm} is used to print summary the result of LSIRM.
#'
#' @param x List; summary of LSIRM with \code{summary.lsirm}.
#' @param ... Additional arguments.
#'
#' @return \code{print.summary.lsirm} return a summary of LSIRM.
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl())
#' summary(lsirm_result)
#' }
#' @rdname print.summary.lsirm
#' @export
print.summary.lsirm <- function(x, ...){
  cat("==========================","\n")
  cat("Summary of model","\n")
  cat("==========================","\n\n")
  cat("Call:\t", deparse(x$call), "\n")
  cat("Model:\t", x$method, "\n")
  cat("Data type:\t", x$dtype, "\n")
  cat("Variable Selection:\t", x$ss, "\n")
  cat("Missing:\t", x$missing, "\n")
  cat(sprintf("MCMC sample of size %i, after burnin of %i iteration",
              x$mcmc.opt$niter,x$mcmc.opt$nburn),"\n\n")
  if(x$n.chains == 1){
    cat("Covariate coefficients posterior means: ","\n\n")
  }else{
    cat("Covariate coefficients posterior means of chain", x$chain,": ","\n\n")
  }
  printCoefmat(x$coef)
  cat("\n---------------------------","\n\n")
  cat("Overall BIC (Smaller is better) :",x$BIC,"\n")
  cat("\nMaximum Log-posterior Iteration: ", "\n")
  printCoefmat(x$map.inf)
}
