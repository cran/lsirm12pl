#' Diagnostic the result of LSIRM model
#'
#' @description \code{diagnostic} is used to diagnostic the result of LSIRM model.
#'
#' @param object object of class \code{lsirm}.
#' @param plot If \code{TRUE}, MCMC diagnostic plots are returned
#' @param draw.item Select items for diagnosis. A default "first" is drawing the first item for selected parameters in \code{which.draw}. Parameter "alpha", however, uses "second" as a default value because the first alpha has an estimation issue. Positions and item names are supported. For instance, if the name of item is "i1", "i2", "i3" and its positions is in order, the result of \code{beta = c("i1","i2","i3")} and \code{beta = c(1,2,3)} are equivalent.
#' @param which.draw Select parameters for diagnosis. For the 1PL model, "beta", "theta" and "gamma" are available. "alpha" is only available in the 2PL model.
#'
#' @return \code{diagnostic} returns plots for checking MCMC convergence for selected parameters.
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5), ncol=10, nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = FALSE, fixed_gamma = FALSE))
#'
#' # 1PL model
#' diagnostic(lsirm_result, plot=TRUE,
#'            which.draw=c("beta","theta","gamma"))
#'
#' # 1PL model, multiple items
#' diagnostic(lsirm_result, plot=TRUE, draw.item=list(beta = c(1,2,3)),
#'            which.draw=c("beta", "gamma"))
#'
#' lsirm_result <- lsirm(data ~ lsirm2pl(spikenslab = FALSE, fixed_gamma = FALSE))
#'
#' # 2PL model
#' diagnostic(lsirm_result, plot=TRUE,
#'            which.draw=c("beta", "theta", "alpha", "gamma"))
#' }

#' @export diagnostic
diagnostic <- function(object, plot=TRUE,
                       draw.item=list(beta="first",
                                      theta="first",
                                      alpha="second"),
                       which.draw=c("beta","gamma")){
  UseMethod("diagnostic")
}

#' @export
diagnostic.lsirm <- function(object, plot=TRUE,
                             draw.item=list(beta="first",
                                            theta="first",
                                            alpha="second"),
                             which.draw=c("beta","gamma"))
{
  mcmc.beta = mcmc.gamma = mcmc.theta = mcmc.alpha = NULL
  if("beta" %in% which.draw) mcmc.beta <- as.mcmc(object$beta)
  if("gamma" %in% which.draw) mcmc.gamma <- as.mcmc(object$gamma)
  if("theta" %in% which.draw) mcmc.theta <- as.mcmc(object$theta)
  if("alpha" %in% which.draw & object$method == "lsirm2pl") mcmc.alpha <- as.mcmc(object$alpha)

  if((is.null(mcmc.beta)) & ("beta" %in% which.draw)) stop("MCMC sample was not stored")
  if((is.null(mcmc.gamma)) & ("gamma" %in% which.draw)) stop("MCMC sample was not stored")
  if((is.null(mcmc.theta)) & ("theta" %in% which.draw)) stop("MCMC sample was not stored")
  if((is.null(mcmc.alpha)) & ("alpha" %in% which.draw)) stop("MCMC sample was not stored. The parameter alpha is only stored in method lsirm2pl")

  if(is.null(draw.item$beta)&("beta" %in% which.draw)) draw.item$beta = "first"
  if(is.null(draw.item$theta)&("theta" %in% which.draw)) draw.item$theta = "first"
  if(is.null(draw.item$alpha)&("alpha" %in% which.draw)) draw.item$alpha = "second"

  acf.val = acf.print = list()
  if((!is.null(mcmc.beta)) & ("beta" %in% which.draw)){
    if(is.null(colnames(object$data))){
      colnames(mcmc.beta) = 1:ncol(mcmc.beta)
    }else{colnames(mcmc.beta) <- colnames(object$data)}
    if(is.numeric(draw.item$beta)){
      draw.item.beta = colnames(mcmc.beta)[draw.item$beta]
    }
    res = list()
    if((length(draw.item$beta) == 1)&(draw.item$beta[1] == "first")){
      draw.item.beta = colnames(mcmc.beta)[1]
    }else{
      draw.item.beta = draw.item$beta
    }
    res$mcmc.sample.beta = mcmc.beta[,draw.item.beta]
    for(i in 1:length(draw.item.beta)){
      ACF = acf(mcmc.beta[,draw.item.beta[i]], plot=F)
      acf.val[[i]] = ACF
      acf.print[[i]] = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
    }
    res$acf.val.beta = acf.val
    res$acf.print.beta = acf.print
    # plot
    if(plot){
      if(length(draw.item.beta) == 1){
        par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
        ts.plot(res$mcmc.sample.beta, ylab="Beta", main = "Trace Plot")
        plot(res$acf.val.beta[[1]], main = "")
        title("ACF", line = 1.6)
        # plot(density(res$mcmc.sample.beta), main = paste0("Density ", draw.item.beta))
        # rug(jitter(res$mcmc.sample.beta))
        plot(res$mcmc.sample.beta, density = T, trace = F, main = "Density", auto.layout = F)
        mtext(paste0("Beta ", draw.item.beta),line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
        par(ask = dev.interactive())
      } else{
        for(i in 1:length(draw.item.beta)){
          if(i > 1){par(ask = dev.interactive())}
          par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
          ts.plot(res$mcmc.sample.beta[,i], ylab="Beta", main = "Trace Plot")
          plot(res$acf.val.beta[[i]],main="")
          title("ACF", line = 1.6)
          # plot(density(res$mcmc.sample.beta[,i]), main = paste0("Density ", draw.item.beta[i]))
          # plot(density(res$mcmc.sample.beta[,i]), main = "Density")
          # rug(jitter(res$mcmc.sample.beta[,i]))
          plot(res$mcmc.sample.beta[,i], density = T, trace = F, main = "Density", auto.layout = F)
          mtext(paste0("Beta ", draw.item.beta[i]),line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
        }
      }
    }
  }

  acf.val = acf.print = list()
  if((!is.null(mcmc.theta)) & ("theta" %in% which.draw)){
    if(is.null(rownames(object$data))){
      colnames(mcmc.theta) = 1:ncol(mcmc.theta)
    }else{colnames(mcmc.theta) <- rownames(object$data)}

    if(is.numeric(draw.item$theta)){
      draw.item.theta = colnames(mcmc.theta)[draw.item$theta]
    }
    res = list()
    if((length(draw.item$theta) == 1)&(draw.item$theta[1] == "first")){
      draw.item.theta = colnames(mcmc.theta)[1]
    }else{
      draw.item.theta = draw.item$theta
    }
    res$mcmc.sample.theta = mcmc.theta[,draw.item.theta]
    for(i in 1:length(draw.item.theta)){
      ACF = acf(mcmc.theta[,draw.item.theta[i]],plot=F)
      acf.val[[i]] = ACF
      acf.print[[i]] = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
    }
    res$acf.val.theta = acf.val
    res$acf.print.theta = acf.print

    # plot
    if(plot){
      if(length(draw.item.theta) == 1){
        par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
        ts.plot(res$mcmc.sample.theta, ylab="Theta", main = "Trace Plot")
        plot(res$acf.val.theta[[1]], main = "")
        title("ACF", line = 1.6)
        # plot(density(res$mcmc.sample.theta), main = paste0("Density ", draw.item.theta))
        # plot(density(res$mcmc.sample.theta), main = "Density")
        plot(res$mcmc.sample.theta, density = T, trace = F, main = "Density", auto.layout = F)
        # rug(jitter(res$mcmc.sample.theta))
        mtext(paste0("Theta ", draw.item.theta),line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
        par(ask = dev.interactive())
      } else{
        for(i in 1:length(draw.item.theta)){
          if(i > 1){par(ask = dev.interactive())}
          par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
          ts.plot(res$mcmc.sample.theta[,i], ylab="Theta", main = "Trace Plot")
          plot(res$acf.val.theta[[i]],main="")
          title("ACF", line = 1.6)
          plot(density(res$mcmc.sample.theta[,i]), main = "Density", auto.layout = F)
          # plot(density(res$mcmc.sample.theta[,i]), main = paste0("Density ", draw.item.theta[i]))
          # rug(jitter(res$mcmc.sample.theta[,i]))
          mtext(paste0("Theta ", draw.item.theta[i]),line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
        }
      }
    }
  }

  acf.val = acf.print = list()
  if((!is.null(mcmc.alpha)) & ("alpha" %in% which.draw)){
    if("first" %in% draw.item$alpha | 1 %in% draw.item$alpha){
      stop("The value of first alpha is fixed to 1 because of identifiability issue.")
    }else{
      if(is.null(colnames(object$data))){
        colnames(mcmc.alpha) = 1:ncol(mcmc.alpha)
      }else{colnames(mcmc.alpha) <- colnames(object$data)}

      if(is.numeric(draw.item$alpha)){
        draw.item.alpha = colnames(mcmc.alpha)[draw.item$alpha]
      }
      res = list()
      if((length(draw.item$alpha) == 1)&(draw.item$alpha[1] == "second")){
        draw.item.alpha = colnames(mcmc.alpha)[2]
      }else{
        draw.item.alpha = draw.item$alpha
      }
      res$mcmc.sample.alpha = mcmc.alpha[,draw.item.alpha]
      for(i in 1:length(draw.item.alpha)){
        ACF = acf(mcmc.alpha[,draw.item.alpha[i]],plot=F)
        acf.val[[i]] = ACF
        acf.print[[i]] = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
      }
      res$acf.val.alpha = acf.val
      res$acf.print.alpha = acf.print

      # plot
      if(plot){
        if(length(draw.item.alpha) == 1){
          par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
          ts.plot(res$mcmc.sample.alpha, ylab="Alpha", main = "Trace Plot")
          plot(res$acf.val.alpha[[1]], main = "")
          title("ACF", line = 1.6)
          # plot(density(res$mcmc.sample.alpha), main = paste0("Density ", draw.item.alpha))
          # plot(density(res$mcmc.sample.alpha), main = "Density")
          # rug(jitter(res$mcmc.sample.alpha))
          plot(res$mcmc.sample.alpha, density = T, trace = F, main = "Density", auto.layout = F)
          mtext(paste0("Alpha ", draw.item.alpha),line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
          par(ask = dev.interactive())
        } else{
          for(i in 1:length(draw.item.alpha)){
            if(i > 1){par(ask = dev.interactive())}
            par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
            ts.plot(res$mcmc.sample.alpha[,i], ylab="Alpha", main = "Trace Plot")
            plot(res$acf.val.alpha[[i]], main = "")
            title("ACF", line = 1.6)
            # plot(density(res$mcmc.sample.alpha[,i]), main = paste0("Density ", draw.item.alpha[i]))
            # plot(density(res$mcmc.sample.alpha[,i]), main = "Density")
            # rug(jitter(res$mcmc.sample.alpha[,i]))
            plot(res$mcmc.sample.alpha[,i], density = T, trace = F, main = "Density", auto.layout = F)
            mtext(paste0("Alpha ", draw.item.alpha[i]),line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
          }
        }
      }
    }
  }

  if((!is.null(mcmc.gamma)) & ("gamma" %in% which.draw)){
    res$mcmc.sample.gamma = mcmc.gamma
    res$acf.val.gamma = acf(res$mcmc.sample.gamma,plot=F)
    res$acf.print.gamma = setNames(drop(res$acf.val.gamma$acf), format(drop(res$acf.val.gamma$lag), digits = 3L))

    # plot
    if(plot){
      par(mfrow=c(1,3), oma=c(0, 0, 4, 0))
      ts.plot(res$mcmc.sample.gamma, ylab="Gamma", main = "Trace Plot")
      plot(res$acf.val.gamma, main = "")
      title("ACF", line = 1.6)
      # plot(density(res$mcmc.sample.gamma), main = "Density")
      # rug(jitter(res$mcmc.sample.gamma))
      plot(res$mcmc.sample.gamma, density = T, trace = F, main = "Density", auto.layout = F)
      mtext("Gamma",line=0.5,side=3,outer=TRUE, cex = 1.5) # add title
      if(("gamma" %in% which.draw)&(length(which.draw)==1)) par(ask = dev.interactive())
    }
  }
  par(mfrow=c(1,1), oma=c(0, 0, 0, 0))
  devAskNewPage(F)
}
