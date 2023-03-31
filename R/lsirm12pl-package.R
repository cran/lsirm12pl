#' lsirm12pl-package
#'
#' Analysis of dichotomous and continuous response data using latent factor by both 1PL LSIRM and 2PL LSIRM as described in Jeon et al. (2021) <doi:10.1007/s11336-021-09762-5>. It includes original 1PL LSIRM and 2PL LSIRM provided for binary response data and its extension for continuous response data. Bayesian model selection with spike-and-slab prior and method for dealing data with missing value under missing at random, missing completely at random are also supported. Various diagnostic plots are available to inspect the latent space and summary of estimated parameters.
#'
#' @docType package
#' @name lsirm12pl
#' @importFrom Rcpp evalCpp
#' @importFrom MCMCpack procrustes
#' @importFrom grDevices boxplot.stats dev.interactive devAskNewPage rainbow
#' @importFrom graphics mtext par rug title
#' @importFrom stats acf density printCoefmat quantile setNames ts.plot
#' @importFrom coda as.mcmc
#' @importFrom spatstat.geom owin
#' @importFrom spatstat.random rpoispp
#' @importFrom plotly ggplotly
#' @importFrom stats median pnorm rnorm runif
#' @importFrom utils capture.output setTxtProgressBar txtProgressBar
#' @import ggplot2 GPArotation dplyr pROC spatstat
#' @useDynLib lsirm12pl
#'
NULL
