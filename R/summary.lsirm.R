#' Summary the result of LSIRM model
#'
#' @description \link{summary} is used to summary the result of LSIRM model.
#'
#' @param object object of class \code{lsirm}.
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
#' \item{dtype}{Type of input data(Binary or Continuous).}
#' \item{ss}{\code{TRUE} if using spike-slab prior}
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
#'
#' @export
summary.lsirm <- function(object, ...)
{
  if(object$method == "lsirm1pl") method = "lpl LSIRM"
  if(object$method == "lsirm2pl") method = "2pl LSIRM"

  res <- list(call = object$call,
              coef = object$beta_summary,
              mcmc.opt = object$mcmc_inf,
              map.inf = object$map_inf,
              BIC = object$bic,
              method = method,
              missing = object$missing,
              dtype = object$dtype,
              ss = object$varselect)
  class(res) <- "summary.lsirm"
  res
}
