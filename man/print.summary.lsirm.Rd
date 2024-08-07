% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.summary.lsirm.R
\name{print.summary.lsirm}
\alias{print.summary.lsirm}
\title{Print the summary the result of LSIRM}
\usage{
\method{print}{summary.lsirm}(x, ...)
}
\arguments{
\item{x}{List; summary of LSIRM with \code{summary.lsirm}.}

\item{...}{Additional arguments.}
}
\value{
\code{print.summary.lsirm} return a summary of LSIRM.
}
\description{
\link{print.summary.lsirm} is used to print summary the result of LSIRM.
}
\examples{
\donttest{
# generate example item response matrix
data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
lsirm_result <- lsirm(data ~ lsirm1pl())
summary(lsirm_result)
}
}
