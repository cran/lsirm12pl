#' Plotting the interaction map of fitted LSIRM model
#'
#' @description \link{plot} is used to plot the latent space of fitted LSIRM model.
#'
#' @param x object of class \code{lsirm1pl}, \code{lsirm2pl}.
#' @param rotation Logical; If TRUE the latent positions are visualized after oblique (oblimin) rotation.
#' @param option character; If value is "interaction", draw the interaction map that represents interactions between respondents and items. If value is "beta", draw the boxplot for the posterior samples of beta. If value is "theta", draw the distribution of the theta estimates per total test score for the \code{data}. If value is "alpha", draw the boxplot for the posterior samples of alpha. The "alpha" is only available for 2pl LSIRM.
#' @param ... Additional arguments for the corresponding function.
#'
#' @return \code{plot} returns the interaction map or boxplot for parameter estimate.
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = FALSE, fixed_gamma = FALSE))
#' plot(lsirm_result)
#'
#' # use oblique rotation
#' plot(lsirm_result, rotation = TRUE)
#' }
#'
#' @export
plot.lsirm = function(x, option = "interaction", rotation=FALSE, ...){
  # x <- NULL
  value <- NULL
  data = x$data
  if(option == "interaction"){
    axis1 <- NULL; axis2 <- NULL
    # if(sum(x$z_estimate) == 0){
    #   item_position = x$ls_mean_item
    #   resp_position = x$ls_mean_respondent
    #   notation = c('v','u')
    # }else{
    #   item_position = x$w_estimate
    #   resp_position = object$z_estimate
    #   notation = c('w','z')
    # }

    item_position = x$w_estimate
    resp_position = x$z_estimate
    notation = c('w','z')

    if((dim(item_position)[2]!=2)|(dim(resp_position)[2]!=2)){
      stop('\"plot_latent\" is implemented for two-dimensional latent space.')
    }

    if(rotation){
      rot <- oblimin(item_position)
      item_position = rot$loadings # item_position <- item_position %*% t(solve(rot$Th))
      resp_position <- resp_position %*% t(solve(rot$Th))
    }

    df1=as.data.frame(item_position)
    df2=as.data.frame(resp_position)
    df1[,3]=notation[1]
    df2[,3]=notation[2]

    colnames(df1)=c('axis1','axis2','source')
    colnames(df2)=c('axis1','axis2','source')

    df=rbind(df2,df1)
    colnames(df)=c('axis1','axis2','source')

    max_coordinate = sapply(df[,c(1,2)], max, na.rm = TRUE)
    min_coordinate = sapply(df[,c(1,2)], min, na.rm = TRUE)
    axis_value = max(abs(c(max_coordinate,min_coordinate)))
    axis_range = c(-axis_value,axis_value)*1.1

    ggplot() +
      geom_point(data = df2, aes(x = axis1, y = axis2), size = 0.7) +
      geom_text(data = df1, aes(x = axis1, y = axis2, label=1:nrow(df1)), color = "red", size = 4, fontface = "bold") +
      xlim(axis_range)+ylim(axis_range) + coord_cartesian(expand = FALSE) + theme_bw() +
      # geom_point(size = 2,alpha=1,shape=16) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"),
            axis.text=element_text(size=16),
            axis.title=element_text(size=14,face="bold"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.title=element_blank(),
            legend.position = c(0.9,0.9),
            legend.text = element_text(size=16),
            plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
      ggtitle("Interaction Map")
  }else if(option == "beta"){
    beta_dataframe <- data.frame(x = rep(1:ncol(x$beta), each= nrow(x$beta)),
                                 value = as.vector(x$beta))
    #outlier
    lower <- min(data.frame(beta_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
    upper <- max(data.frame(beta_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5, 2])
    ylim1 = c(lower, upper)

    if(ncol(data) > 30){
      ggplot(data=beta_dataframe, aes(x=x,y=value, group=x)) + geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = ylim1*1.05) +
        scale_x_continuous(breaks = round(seq(from = 0, to = ncol(x$beta), length.out = 10))) +
        xlab("Item number") + ylab("Beta estimates")+
        theme(axis.text.x = element_text(face="bold",size=13),
              axis.text.y = element_text( face="bold",size=15),
              axis.title=element_text(size=15, face='bold'))
    }else{
      ggplot(data=beta_dataframe, aes(x=x,y=value, group=x))+ geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = ylim1*1.05) +
        scale_x_continuous(breaks = 0:ncol(x$beta)) +
        xlab("Item number") + ylab("Beta estimates")+
        theme(axis.text.x = element_text(face="bold",size=13),
              axis.text.y = element_text( face="bold",size=15),
              axis.title=element_text(size=15, face='bold'))
    }
  }else if(option == "theta"){
    if(exists(x$missing.val)){
    data[data==x$missing.val] = NA
    }

    total_score <- rowSums(data, na.rm = TRUE)

    if(sum(data != 1 & data != 0, na.rm = TRUE) > 1){
      # Continuous Data
      binning <- cut(total_score, breaks = seq(from = min(total_score), to = max(total_score),
                                               length.out = 11), include.lowest = TRUE)

      theta_temp <- data.frame(x = binning, value = x$theta_estimate)

      # # if null level exist
      # null_level <- levels(binning)[table(binning) == 0]
      # if(length(null_level) != 0) theta_temp <- rbind(theta_temp, data.frame(x = null_level, value = NA))

      ggplot(data=theta_temp, aes(x=x,y=value, group=x))+ geom_boxplot() +
        xlab("Sum score") + ylab("Theta estimates")+
        theme(axis.text.x = element_text(face="bold",size=8),
              axis.text.y = element_text( face="bold",size=15),
              axis.title=element_text(size=15, face='bold'))
    }else{
      # Binary Data
      theta_temp <- data.frame(x = total_score, value = x$theta_estimate)
      if(ncol(data) > 30){
        ggplot(data=theta_temp, aes(x=x,y=value, group=x))+ geom_boxplot() +
          xlab("Sum score") + ylab("Theta estimates")+
          scale_x_continuous(breaks = round(seq(from = 0, to = ncol(x$beta), length.out = 10))) +
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold'))
      }else{
        ggplot(data=theta_temp, aes(x=x, y=value, group = x))+ geom_boxplot() +
          xlab("Sum score") + ylab("Theta estimates")+
          scale_x_continuous(breaks = 0:ncol(x$beta)) +
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold'))
      }
    }
  }else if(option == "alpha"){
    if(is.null(x$alpha) == TRUE){
      stop('The option "alpha" is only available for 2pl LSIRM.')
    }else{
      alpha_dataframe <- data.frame(x = rep(1:ncol(x$alpha), each= nrow(x$alpha)),
                                    value = as.vector(x$alpha))

      #outlier
      lower <- min(data.frame(alpha_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
      upper <- max(data.frame(alpha_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = 'drop'))[1:ncol(data)*5, 2])
      ylim1 = c(lower, upper)

      if(ncol(data) > 30){
        ggplot(data=alpha_dataframe, aes(x=x,y=value, group=x)) + geom_boxplot(outlier.shape = NA) +
          coord_cartesian(ylim = ylim1*1.05) +
          scale_x_continuous(breaks = round(seq(from = 0, to = ncol(x$alpha), length.out = 10))) +
          xlab("Item number") + ylab("Alpha estimates")+
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold'))
      }else{
        ggplot(data=alpha_dataframe, aes(x=x,y=value, group=x))+ geom_boxplot(outlier.shape = NA) +
          coord_cartesian(ylim = ylim1*1.05) +
          scale_x_continuous(breaks = 0:ncol(x$alpha)) +
          xlab("Item number") + ylab("Alpha estimates")+
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold'))
      }
    }
  }
}


