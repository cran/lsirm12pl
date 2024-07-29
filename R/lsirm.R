#' Fit a LSIRM ( Latent Space Item Response Model)
#'
#' @description
#' \link{lsirm} is used to fit 1PL LSIRM and 2PL LSIRM using Bayesian method as described in Jeon et al. (2021).
#'
#' @param formula The form of formula is \code{lsirm(A ~ <term 1>(<term 2>, <term 3> ...))}, where \code{A} is binary or continuous item response matrix to be analyzed, \code{<term1>} is the model you want to fit and has one of the following values: "lsirm1pl" and "lsirm2pl"., and \code{<term 2>}, \code{<term 3>}, etc. are each option for the model.
#' @param ... Additional arguments for the corresponding function.
#'
#' @details
#' The descriptions of options for each model, such as \code{<term 2>} and \code{<term 3>}, are included in \code{\link{lsirm1pl}} for 1PL LSIRM and \code{\link{lsirm2pl}} for 2PL LSIRM.
#'
#' @return \code{lsirm} returns an object of class \code{list}.
#'
#' See corresponding functions such as \code{\link{lsirm1pl}} for 1PL LSIRM and \code{\link{lsirm2pl}} for 2PL LSIRM.
#'
#' @seealso
#' \code{\link{lsirm1pl}} for 1PL LSIRM.
#'
#' \code{\link{lsirm2pl}} for 2PL LSIRM.
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#'
#' lsirm_result <- lsirm(data~lsirm1pl())
#' lsirm_result <- lsirm(data~lsirm2pl())
#'
#' }
#' @export
lsirm = function(formula, ...) UseMethod("lsirm")



#' Formula function for LSIRM
#'
#' @description \link{lsirm.formula} is formula object.
#'
#' @param formula The form of formula is \code{lsirm(A ~ <term 1>(<term 2>, <term 3> ...))}, where \code{A} is binary or continuous item response matrix to be analyzed, \code{<term1>} is the model you want to fit and has one of the following values: "lsirm1pl" and "lsirm2pl"., and \code{<term 2>}, \code{<term 3>}, etc., are each option for the model.
#' @param ... Additional arguments for the corresponding function.
#'
#' @export
lsirm.formula = function(formula, ...){
  data <- get(eval(formula)[[2]])
  argument <- rlang::call_args(eval(formula)[[3]])
  argument$data <- data
  func <- as.character(eval(formula)[[3]][[1]])
  output <- do.call(func,argument)
  if((!is.null(argument$chains))&(output$chains > 1)){
    for(i in 1:argument$chains){
      output[[i]]$call <- match.call()
      class(output[[i]]) <- "lsirm"
    }
  }else{
    output$call <- match.call()
    class(output) <- "lsirm"
  }
  return(output)
}





#' Fit a 1PL LSIRM for binary and continuous item response data
#' @description
#' \code{\link{lsirm1pl}} integrates all functions related to 1PL LSIRM. Various 1PL LSIRM function can be used by setting the \code{spikenslab}, \code{fixed_gamma}, and \code{missing_data} arguments.
#'
#' This function can be used regardless of the data type, providing a unified approach to model fitting.
#'
#'
#' @param data Matrix; a binary or continuous item response matrix for analysis. Each row represents a respondent, and each column contains responses to the corresponding item.
#' @param spikenslab Logical; specifies whether to use a model selection approach. Default is FALSE.
#' @param fixed_gamma Logical; indicates whether to fix gamma at 1. Default is FALSE.
#' @param missing_data Character; the type of missing data assumed. Options are NA, "mar", or "mcar". Default is NA.
#' @param chains Integer; the number of MCMC chains to run. Default is 1.
#' @param multicore Integer; the number of cores to use for parallel execution. Default is 1.
#' @param seed Integer; the seed number for MCMC fitting. Default is NA.
#' @param ndim Integer; the dimension of the latent space. Default is 2.
#' @param niter Integer; the total number of MCMC iterations to run. Default is 15000.
#' @param nburn Integer; the number of initial MCMC iterations to discard as burn-in. Default is 2500.
#' @param nthin Integer; the number of MCMC iterations to thin. Default is 5.
#' @param nprint Integer; the interval at which MCMC samples are displayed during execution. Default is 500.
#' @param jump_beta Numeric; the jumping rule for the beta proposal density. Default is 0.4.
#' @param jump_theta Numeric; the jumping rule for the theta proposal density. Default is 1.0.
#' @param jump_z Numeric; the jumping rule for the z proposal density. Default is 0.5.
#' @param jump_w Numeric; the jumping rule for the w proposal density. Default is 0.5.
#' @param pr_mean_beta Numeric; the mean of the normal prior for beta. Default is 0.
#' @param pr_sd_beta Numeric; the standard deviation of the normal prior for beta. Default is 1.0.
#' @param pr_mean_theta Numeric; the mean of the normal prior for theta. Default is 0.
#' @param pr_a_theta Numeric; the shape parameter of the inverse gamma prior for the variance of theta. Default is 0.001.
#' @param pr_b_theta Numeric; the scale parameter of the inverse gamma prior for the variance of theta. Default is 0.001.
#' @param \dots Additional arguments for the for various settings. Refer to the functions in the Details.
#'
#' @return \code{lsirm1pl} returns an object of list.
#' The basic return list containing the following components:
#' \item{data}{A data frame or matrix containing the variables used in the model.}
#' \item{bic}{A numeric value representing the Bayesian Information Criterion (BIC).}
#' \item{mcmc_inf}{Details about the number of MCMC iterations, burn-in periods, and thinning intervals.}
#' \item{map_inf}{The log maximum a posteriori (MAP) value and the iteration number at which this MAP value occurs.}
#' \item{beta_estimate}{Posterior estimates of the beta parameter.}
#' \item{theta_estimate}{Posterior estimates of the theta parameter.}
#' \item{sigma_theta_estimate}{Posterior estimates of the standard deviation of theta.}
#' \item{z_estimate}{Posterior estimates of the z parameter.}
#' \item{w_estimate}{Posterior estimates of the w parameter.}
#' \item{beta}{Posterior samples of the beta parameter.}
#' \item{theta}{Posterior samples of the theta parameter.}
#' \item{theta_sd}{Posterior samples of the standard deviation of theta.}
#' \item{z}{Posterior samples of the z parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{w}{Posterior samples of the w parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{accept_beta}{Acceptance ratio for the beta parameter.}
#' \item{accept_theta}{Acceptance ratio for the theta parameter.}
#' \item{accept_z}{Acceptance ratio for the z parameter.}
#' \item{accept_w}{Acceptance ratio for the w parameter.}
#' \item{...}{Additional return values for various settings. Refer to the functions in the Details.}
#'
#' @details
#'  Additional arguments and return values for each function are documented in the respective function's description.
#'
#' * For LSIRM with data included missing value are detailed in \link{lsirm1pl_mar} and \link{lsirm1pl_mcar}.
#'
#' * For LSIRM using the spike-and-slab model selection approach are detailed in \link{lsirm1pl_ss}.
#'
#' * For continuous version of LSIRM are detailed in \link{lsirm1pl_normal_o}.
#'
#' @note If both \code{spikenslab} and \code{fixed_gamma} are set \code{TRUE}, it returns error because both are related to \code{gamma}.
#'
#' @seealso
#' The LSIRM for 1PL LSIRM for binary item response data as following:
#'
#' \code{\link{lsirm1pl_o}}, \code{\link{lsirm1pl_fixed_gamma}}, \code{\link{lsirm1pl_mar}},\code{\link{lsirm1pl_mcar}}, \code{\link{lsirm1pl_fixed_gamma_mar}}, \code{\link{lsirm1pl_fixed_gamma_mcar}}, \code{\link{lsirm1pl_ss}}, \code{\link{lsirm1pl_mar_ss}}, and \code{\link{lsirm1pl_mcar_ss}}
#'
#' The LSIRM for 1PL LSIRM for continuous item response data as following:
#'
#' \code{\link{lsirm1pl_normal_o}}, \code{\link{lsirm1pl_normal_fixed_gamma}}, \code{\link{lsirm1pl_normal_mar}},  \code{\link{lsirm1pl_normal_mcar}},\code{\link{lsirm1pl_normal_fixed_gamma_mar}}, \code{\link{lsirm1pl_normal_fixed_gamma_mcar}}, \code{\link{lsirm1pl_normal_ss}}, \code{\link{lsirm1pl_normal_mar_ss}}, \code{\link{lsirm1pl_normal_mcar_ss}}
#'
#' @details
#' For 1PL LSIRM with binary item response data, the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term: \deqn{logit(P(Y_{j,i} = 1|\theta_j,\beta_i,\gamma,z_j,w_i))=\theta_j+\beta_i-\gamma||z_j-w_i||}
#'
#' For 1PL LSIRM with continuous item response data, the continuous value of response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term: \deqn{Y_{j,i} = \theta_j+\beta_i-\gamma||z_j-w_i|| + e_{j,i}} where the error \eqn{e_{j,i} \sim N(0,\sigma^2)}.
#'
#'@examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm1pl(data)
#'
#'# The code following can achieve the same result.
#' lsirm_result <- lsirm(data~lsirm1pl())
#'
#' }
#' @export
lsirm1pl = function(data, spikenslab = FALSE, fixed_gamma = FALSE, missing_data = NA, chains = 1, multicore = 1, seed = NA,
                    ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta, jump_z, jump_w, pr_mean_beta,
                    pr_sd_beta, pr_mean_theta, pr_a_theta, pr_b_theta, ...) {
  if(!is.na(seed)){
    set.seed(seed)
  }
  cat("\n Fitting LSIRM with MCMC algorithm\n")
  if(spikenslab == FALSE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(sum(is.na(data)) > 0){
      stop("The data contains NA. Set missing_data to 'mcar' or 'mar'.")
    }
    if(check.datatype(data, ...)){ # if missing.val imputed
      ####lsirm1pl_o -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_o(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_o(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    } else{
      ####lsirm1pl_normal_o -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_o(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_o(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_o(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      ####lsirm1pl_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_mar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("\n\n Chain %d / %d completed \n\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"

          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_mar(..., data)
        output$dtype <- "binary"

      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm1pl_normal_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_mar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_mar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      ####lsirm1pl_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_mcar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()


          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_mcar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    } else{
      ####lsirm1pl_normal_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_mcar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_mcar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      ####lsirm1pl_fixed_gamma -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_fixed_gamma(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm1pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_fixed_gamma(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_fixed_gamma(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    } else{
      ####lsirm1pl_normal_fixed_gamma -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_fixed_gamma(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_fixed_gamma(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_fixed_gamma(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      ####lsirm1pl_fixed_gamma_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_fixed_gamma_mar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm1pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_fixed_gamma_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_fixed_gamma_mar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm1pl_normal_fixed_gamma_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_fixed_gamma_mar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_fixed_gamma_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_fixed_gamma_mar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      ####lsirm1pl_fixed_gamma_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_fixed_gamma_mcar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm1pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_fixed_gamma_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_fixed_gamma_mcar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm1pl_normal_fixed_gamma_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_fixed_gamma_mcar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_fixed_gamma_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_fixed_gamma_mcar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      ####lsirm1pl_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_ss(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm1pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_ss(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm1pl_normal_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_ss(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_ss(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      ####lsirm1pl_mar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_mar_ss(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm1pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_mar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_mar_ss(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm1pl_normal_mar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_mar_ss(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_mar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_mar_ss(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      ####lsirm1pl_mcar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_mcar_ss(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm1pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_mcar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_mcar_ss(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm1pl_normal_mcar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm1pl_normal_mcar_ss(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm1pl_normal_mcar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm1pl_normal_mcar_ss(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else{
    stop('The options "spikenslab" and "fixed_gamma" cannot be set TRUE at the same time.')
  }
  if(chains > 1){
    for(i in 1:chains){
      output[[i]]$call <- match.call()
      output$method <- "lsirm1pl"
      output[[i]]$method <- "lsirm1pl"
      output[[i]]$missing <- missing_data
      output[[i]]$varselect <- spikenslab
      output$chains = chains
      class(output) <- "lsirm"

    }
  }else{
    output$call <- match.call()
    output$method <- "lsirm1pl"
    output$missing <- missing_data
    output$varselect <- spikenslab
    output$chains = chains
    class(output) <- "lsirm"
  }
  return(output)
}



#'  Fit a 2pl LSIRM for binary and continuous item resopnse data
#'
#' @description
#' \code{\link{lsirm2pl}} integrates all functions related to 2PL LSIRM. Various 2PL LSIRM function can be used by setting the \code{spikenslab}, \code{fixed_gamma}, and \code{missing_data} arguments.
#'
#' This function can be used regardless of the data type, providing a unified approach to model fitting.
#'
#' @param data Matrix; a binary or continuous item response matrix for analysis. Each row represents a respondent, and each column contains responses to the corresponding item.
#' @param spikenslab Logical; specifies whether to use a model selection approach. Default is FALSE.
#' @param fixed_gamma Logical; indicates whether to fix gamma at 1. Default is FALSE.
#' @param missing_data Character; the type of missing data assumed. Options are NA, "mar", or "mcar". Default is NA.
#' @param chains Integer; the number of MCMC chains to run. Default is 1.
#' @param multicore Integer; the number of cores to use for parallel execution. Default is 1.
#' @param seed Integer; the seed number for MCMC fitting. Default is NA.
#' @param ndim Integer; the dimension of the latent space. Default is 2.
#' @param niter Integer; the total number of MCMC iterations to run. Default is 15000.
#' @param nburn Integer; the number of initial MCMC iterations to discard as burn-in. Default is 2500.
#' @param nthin Integer; the number of MCMC iterations to thin. Default is 5.
#' @param nprint Integer; the interval at which MCMC samples are displayed during execution. Default is 500.
#' @param jump_beta Numeric; the jumping rule for the beta proposal density. Default is 0.4.
#' @param jump_theta Numeric; the jumping rule for the theta proposal density. Default is 1.0.
#' @param jump_alpha Numeric; the jumping rule for the alpha proposal density. Default is 1.0.
#' @param jump_z Numeric; the jumping rule for the z proposal density. Default is 0.5.
#' @param jump_w Numeric; the jumping rule for the w proposal density. Default is 0.5.
#' @param pr_mean_beta Numeric; the mean of the normal prior for beta. Default is 0.
#' @param pr_sd_beta Numeric; the standard deviation of the normal prior for beta. Default is 1.0.
#' @param pr_mean_theta Numeric; the mean of the normal prior for theta. Default is 0.
#' @param pr_mean_alpha Numeric; the mean of the log normal prior for alpha. Default is 0.5.
#' @param pr_sd_alpha Numeric; the standard deviation of the log normal prior for alpha. Default is 1.0.
#' @param pr_a_theta Numeric; the shape parameter of the inverse gamma prior for the variance of theta. Default is 0.001.
#' @param pr_b_theta Numeric; the scale parameter of the inverse gamma prior for the variance of theta. Default is 0.001.
#' @param \dots Additional arguments for the for various settings. Refer to the functions in the Details.
#'
#' @return \code{lsirm2pl} returns an object of list.
#' The basic return list containing the following components:
#' \item{data}{A data frame or matrix containing the variables used in the model.}
#' \item{bic}{A numeric value representing the Bayesian Information Criterion (BIC).}
#' \item{mcmc_inf}{Details about the number of MCMC iterations, burn-in periods, and thinning intervals.}
#' \item{map_inf}{The log maximum a posteriori (MAP) value and the iteration number at which this MAP value occurs.}
#' \item{beta_estimate}{Posterior estimates of the beta parameter.}
#' \item{theta_estimate}{Posterior estimates of the theta parameter.}
#' \item{sigma_theta_estimate}{Posterior estimates of the standard deviation of theta.}
#' \item{alpha_estimate}{posterior estimates of alpha parameter..}
#' \item{z_estimate}{Posterior estimates of the z parameter.}
#' \item{w_estimate}{Posterior estimates of the w parameter.}
#' \item{beta}{Posterior samples of the beta parameter.}
#' \item{theta}{Posterior samples of the theta parameter.}
#' \item{theta_sd}{Posterior samples of the standard deviation of theta.}
#' \item{alpha}{Posterior samples of the alpha parameter.}
#' \item{z}{Posterior samples of the z parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{w}{Posterior samples of the w parameter, represented as a 3-dimensional matrix where the last axis denotes the dimension of the latent space.}
#' \item{accept_beta}{Acceptance ratio for the beta parameter.}
#' \item{accept_theta}{Acceptance ratio for the theta parameter.}
#' \item{accept_z}{Acceptance ratio for the z parameter.}
#' \item{accept_w}{Acceptance ratio for the w parameter.}
#' \item{accept_alpha}{Acceptance ratio for the alpha parameter.}
#' \item{...}{Additional return values for various settings. Refer to the functions in the Details.}
#'
#' @details
#' Additional arguments and return values for each function are documented in the respective function's description.
#'
#' * For 2PL LSIRM with data included missing value are detailed in \link{lsirm2pl_mar} and \link{lsirm2pl_mcar}.
#'
#' * For 2PL LSIRM using the spike-and-slab model selection approach are detailed in \link{lsirm2pl_ss}.
#'
#' * For continuous version of 2PL LSIRM are detailed in \link{lsirm2pl_normal_o}.
#'
#' @note If both \code{spikenslab} and \code{fixed_gamma} are set \code{TRUE}, it returns error because both are related to \code{gamma}.
#'
#' @seealso
#' The 2PL LSIRM for binary item response data as following:
#'
#' \code{\link{lsirm2pl_o}}, \code{\link{lsirm2pl_fixed_gamma}}, \code{\link{lsirm2pl_mar}},\code{\link{lsirm2pl_mcar}}, \code{\link{lsirm2pl_fixed_gamma_mar}}, \code{\link{lsirm2pl_fixed_gamma_mcar}}, \code{\link{lsirm2pl_ss}}, \code{\link{lsirm2pl_mar_ss}}, and \code{\link{lsirm2pl_mcar_ss}}
#'
#' The 2PL LSIRM for continuous item response data as following:
#'
#' \code{\link{lsirm2pl_normal_o}}, \code{\link{lsirm2pl_normal_fixed_gamma}}, \code{\link{lsirm2pl_normal_mar}},  \code{\link{lsirm2pl_normal_mcar}},\code{\link{lsirm1pl_normal_fixed_gamma_mar}}, \code{\link{lsirm2pl_normal_fixed_gamma_mcar}}, \code{\link{lsirm2pl_normal_ss}}, \code{\link{lsirm2pl_normal_mar_ss}}, \code{\link{lsirm2pl_normal_mcar_ss}}
#'
#' @details
#' For 2PL LSIRM with binary item response data, the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term. For 2pl model, the the item effect is assumed to have additional discrimination parameter \eqn{\alpha_i} multiplied by \eqn{\theta_j}: \deqn{logit(P(Y_{j,i} = 1|\theta_j,\alpha_i,\beta_i,\gamma,z_j,w_i))=\theta_j*\alpha_i+\beta_i-\gamma||z_j-w_i||}
#'
#' For 2PL LSIRM with continuous item response data, the continuous value of response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term. For 2pl model, the the item effect is assumed to have additional discrimination parameter \eqn{\alpha_i} multiplied by \eqn{\theta_j}: \deqn{Y_{j,i} = \theta_j+\beta_i-\gamma||z_j-w_i|| + e_{j,i}} where the error \eqn{e_{j,i} \sim N(0,\sigma^2)}
#'@examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm2pl(data)
#'
#'# The code following can achieve the same result.
#' lsirm_result <- lsirm(data~lsirm2pl())
#'
#' }
#' @export
lsirm2pl = function(data, spikenslab = FALSE, fixed_gamma = FALSE, missing_data = NA, chains = 1, multicore = 1, seed = NA,
                    ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta, jump_alpha, jump_z, jump_w, pr_mean_beta,
                    pr_sd_beta, pr_mean_theta, pr_a_theta, pr_b_theta,
                    pr_mean_alpha, pr_sd_alpha, ...) {
  if(!is.na(seed)){
    set.seed(seed)
  }
  cat("\n Fitting LSIRM with MCMC algorithm\n")
  if(spikenslab == FALSE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){ # if missing.val imputed
      ####lsirm2pl_o -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_o(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_o(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    } else{
      ####lsirm2pl_normal_o -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_o(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_o(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_o(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      ####lsirm2pl_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_mar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_mar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm2pl_normal_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_mar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_mar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      ####lsirm2pl_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_mcar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()


          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_mcar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    } else{
      ####lsirm2pl_normal_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_mcar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_mcar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      ####lsirm2pl_fixed_gamma -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_fixed_gamma(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_fixed_gamma(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_fixed_gamma(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }
    } else{
      ####lsirm2pl_normal_fixed_gamma -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_fixed_gamma(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_fixed_gamma(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_fixed_gamma(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      ####lsirm2pl_fixed_gamma_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_fixed_gamma_mar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_fixed_gamma_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_fixed_gamma_mar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm2pl_normal_fixed_gamma_mar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_fixed_gamma_mar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_fixed_gamma_mar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_fixed_gamma_mar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      ####lsirm2pl_fixed_gamma_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_fixed_gamma_mcar(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_fixed_gamma_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_fixed_gamma_mcar(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm2pl_normal_fixed_gamma_mcar -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_fixed_gamma_mcar(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_fixed_gamma_mcar(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_fixed_gamma_mcar(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      ####lsirm2pl_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_ss(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_ss(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm2pl_normal_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_ss(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_ss(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      ####lsirm2pl_mar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_mar_ss(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_mar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_mar_ss(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm2pl_normal_mar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_mar_ss(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_mar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_mar_ss(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      ####lsirm2pl_mcar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_mcar_ss(..., data)
            output[[i]]$dtype <- "binary"
            cat(sprintf("Chain %d / %d completed\n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          # clusterExport(cl,c("lsirm2pl_o","chunks"))

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_mcar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            # output <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_o(..., data)},
            #                               data = data, ..., simplify = F),
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )

          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "binary"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_mcar_ss(..., data)
        output$dtype <- "binary"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    } else{
      ####lsirm2pl_normal_mcar_ss -----------
      if(chains > 1){
        if(multicore <= 1){
          output = list()
          for(i in 1:chains){
            output[[i]] <- lsirm2pl_normal_mcar_ss(..., data)
            output[[i]]$dtype <- "continuous"
            cat(sprintf("Chain %d / %d completed \n", i, chains))
          }
        }else{
          if(chains < multicore){
            stop("Error: The number of chains must not be less than the number of cores. Please adjust the number of chains or cores to optimize parallel processing.")
          }else{
            cl <- makeCluster(multicore)
          }

          q <- chains %/% multicore
          r <- chains %% multicore
          chunks <- list()

          for(i in 1:q){
            chunks[[i]] = seq(1,multicore,1)
          }

          if(r != 0){ chunks[[q+1]] = seq(1,r,1) }
          output <- list()

          tryCatch(
            for(c in 1:length(chunks)){
              X <- chunks[[c]]
              output[[c]] <- parallel::parSapply(cl, X, function(X,data,...){lsirm2pl_normal_mcar_ss(..., data)},
                                                 data = data, ..., simplify = F)
              cat(sprintf("\n Chunk %d / %d processed \n", c, length(chunks)))
              # cat(sprintf("\n Core %d / %d finished",c,length(chunks)),"\n")
            },
            error = function(e){
              parallel::stopCluster(cl)
              stop("An error occurred during parallel execution.")
            }
          )
          output <- do.call("c",output)

          for(i in 1:chains){
            output[[i]]$dtype <- "continuous"
          }
          parallel::stopCluster(cl)
        }
      }else if(chains == 1){
        if(chains < multicore){
          warning("Warning: The number of chains is equal to 1. Please adjust the number of chains or cores to optimize parallel processing.")
        }
        output <- lsirm2pl_normal_mcar_ss(..., data)
        output$dtype <- "continuous"
      }else{
        stop("The number of chains must be an integer greater than 1.")
      }

    }

  }else{
    stop('The options "spikenslab" and "fixed_gamma" cannot be set TRUE at the same time.')
  }
  if(chains > 1){
    for(i in 1:chains){
      output[[i]]$call <- match.call()
      output$method <- "lsirm2pl"
      output[[i]]$method <- "lsirm2pl"
      output[[i]]$missing <- missing_data
      output[[i]]$varselect <- spikenslab
      output$chains = chains
      class(output) <- "lsirm"
    }
  }else{
    output$call <- match.call()
    output$method <- "lsirm2pl"
    output$missing <- missing_data
    output$varselect <- spikenslab
    class(output) <- "lsirm"
    output$chains = 1
  }
  return(output)
}




