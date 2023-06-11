#' Goodness-of-fit LSIRM model
#'
#' @description \link{gof} is goodness-of-fit the latent space of fitted LSIRM model.
#'
#' @param object object of class \code{lsirm}.
#' @param nsim Integer; Number of simulation. Default is 500
#'
#' @return \code{gof} returns the boxplot or AUC plot
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = FALSE, fixed_gamma = FALSE))
#' gof(lsirm_result)
#' }

#' @export gof
gof <- function(object, nsim=500){
  UseMethod("gof")
}

#' @export
gof.lsirm = function(object, nsim=500){
  ind <- NULL
  values <- NULL
  x <- NULL
  y <- NULL

  data = object$data
  e.sigma = ifelse(object$dtype=="continuous",object$sigma_estimate,1)
  gamma = ifelse(!is.null(object$gamma_estimate),object$gamma_estimate,1)
  P = make_prob(nsim, beta=object$beta_estimate,theta=object$theta_estimate,
                gamma=gamma, sigma=e.sigma,
                item = object$w_estimate, res = object$z_estimate,
                type = object$dtype)

  colnames(P) = c(1:ncol(data))
  dat = stack(as.data.frame(P))
  colnames(dat) <- c("values", "ind")

  data_temp <- data.frame(x = c(1:ncol(data)), y = colMeans(data))

  if(ncol(data)>30){
    gofp = ggplot(dat) +
      geom_boxplot(aes(x = ind, y = values), outlier.shape = NA)+
      scale_x_discrete(breaks = round(seq(from = 0, to = ncol(data), length.out = 10))) +
      geom_point(data = data_temp,
                 aes(x=x, y=y),
                 color="red",
                 size=1)+
      xlab("Item number")+
      theme(axis.text.x = element_text(face="bold",size=13),
            axis.text.y = element_text( face="bold",size=15),
            axis.title = element_text(size=15, face='bold'),
            plot.margin = margin(1,1,1.5,1.2,"cm"))
  }else{
    gofp = ggplot(dat) +
      geom_boxplot(aes(x = ind, y = values), outlier.shape = NA)+
      geom_point(data = data_temp,
                 aes(x=x, y=y),
                 color="red",
                 size=1)+
      xlab("Item number")+
      theme(axis.text.x = element_text(face="bold",size=13),
            axis.text.y = element_text( face="bold",size=15),
            axis.title = element_text(size=15, face='bold'),
            plot.margin = margin(1,1,1.5,1.2,"cm"))
  }
  if(object$dtype == "continuous"){
    gofp = gofp+
      labs(title="Goodness of fit")+
      theme(title = element_text(size=25,face="bold",hjust=0.5),
            plot.title = element_text(hjust = 0.5))
    return(gofp)
  }else{
    if(object$method == "lsirm1pl"){
      aucp = roc_1pl(object)
      title <- grid::textGrob("Goodness of fit",
                              gp=gpar(font=2,fontsize=25))

      # Add a zeroGrob of height 2cm on top of the title
      title <- arrangeGrob(zeroGrob(), title,
                           widths = unit(1, 'npc'),
                           heights = unit(c(1, 1), c('cm', 'npc')),
                           as.table = FALSE)
      gridExtra::grid.arrange(gofp, aucp, ncol=2,
                              top=title)
    }else{
      aucp = roc_2pl(object)
      title <- grid::textGrob("Goodness of fit",
                              gp=gpar(font=2,fontsize=25))

      # Add a zeroGrob of height 2cm on top of the title
      title <- arrangeGrob(zeroGrob(), title,
                           widths = unit(1, 'npc'),
                           heights = unit(c(1, 1), c('cm', 'npc')),
                           as.table = FALSE)
      gridExtra::grid.arrange(gofp, aucp, ncol=2,
                              top=title)
    }
  }
}

make_prob = function(nsim,beta,theta,gamma,sigma=1,item,res,type="continuous"){
  cat("\nSimulation Start\n\n")
  pb <- txtProgressBar(title = "progress bar", min = 0, max = nsim,
                       style = 3, width = 50)
  n.i = nrow(item)
  n.r = nrow(res)
  S = matrix(nrow=nsim,ncol=n.i)
  P = matrix(nrow=n.r,ncol=n.i)
  if(type=="binary"){
    for(i in 1:n.i){
      for(j in 1:n.r){
        d = sum((item[i,] - res[j,])^2)^{1/2}
        P[j,i] = exp(beta[i]+theta[j]-gamma*d)/(1+exp(beta[i]+theta[j]-gamma*d))
      }
    }
    for(k in 1:nsim){
      setTxtProgressBar(pb, k, label = paste( round(k/nsim * 100, 0), "% done"))
      S[k,] = colMeans(matrix(rbinom(n.r*n.i,1,P),nrow=n.r,ncol=n.i))
    }
  }else if(type=="continuous"){
    for(i in 1:n.i){
      for(j in 1:n.r){
        d = sum((item[i,] - res[j,])^2)^{1/2}
        P[j,i] = beta[i]+theta[j]-gamma*d
      }
    }
    for(k in 1:nsim){
      setTxtProgressBar(pb, k, label = paste( round(k/nsim * 100, 0), "% done"))
      S[k,] = colMeans(P + matrix(rnorm(n.i*n.r,0,sigma),nrow=n.r,ncol=n.i))
    }
  }
  return(S)
}
