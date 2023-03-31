#' Summary the result of LSIRM model
#'
#' @description \link{summary.lsirm} is used to summary the result of LSIRM model.
#'
#' @param x object of class \code{lsirm1pl}, \code{lsirm2pl}.
#' @param \dots Additional arguments.
#'
#' @return \code{summary.lsirm} contains following elements. A print method is available.
#' \item{call}{R call used to fit the model.}
#' \item{coef}{Covariate coefficients posterior means.}
#' \item{mcmc.opt}{The number of mcmc iteration, burn-in periods, and thinning intervals.}
#' \item{map.inf}{value of log maximum a posterior and iteration number which have log maximum a posterior.}
#' \item{BIC}{Numeric value with the corresponding BIC.}
#' \item{method}{1pl LSIRM or 2pl LSIRM }
#' \item{missing}{The assumed missing type. One of NA, "mar" and "mcar". Default uses NA.}
#' \item{dtype}{Binary or Continuous}
#' \item{ss}{\code{TRUE} if using spike-slab prior}
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl())
#' summary(lsirm_result)
#' }
#'
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
  cat("Covariate coefficients posterior means: ","\n\n")
  printCoefmat(x$coef)
  cat("\n---------------------------","\n\n")
  cat("Overall BIC (Smaller is better) :",x$BIC,"\n")
  cat("\nMaximum Log-posterior Iteration: ", "\n")
  printCoefmat(x$map.inf)
}
