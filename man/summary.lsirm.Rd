% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.lsirm.R
\name{summary.lsirm}
\alias{summary.lsirm}
\title{Summary the result of LSIRM model}
\usage{
\method{summary}{lsirm}(object, ...)
}
\arguments{
\item{object}{object of class \code{lsirm}.}

\item{\dots}{Additional arguments.}
}
\value{
\code{summary.lsirm} contains following elements. A print method is available.
\item{call}{R call used to fit the model.}
\item{coef}{Covariate coefficients posterior means.}
\item{mcmc.opt}{The number of mcmc iteration, burn-in periods, and thinning intervals.}
\item{map.inf}{value of log maximum a posterior and iteration number which have log maximum a posterior.}
\item{BIC}{Numeric value with the corresponding BIC.}
\item{method}{1pl LSIRM or 2pl LSIRM }
\item{missing}{The assumed missing type. One of NA, "mar" and "mcar". Default uses NA.}
\item{dtype}{Type of input data(Binary or Continuous).}
\item{ss}{\code{TRUE} if using spike-slab prior}
}
\description{
\link{summary} is used to summary the result of LSIRM model.
}
\examples{
\donttest{
# generate example item response matrix
data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)

# 1PL LSIRM object
lsirm_result <- lsirm(data ~ lsirm1pl())
summary(lsirm_result)
}

}
