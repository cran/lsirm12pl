#' boxplot of parameters of fitted LSIRM model
#' 
#' @description \link{plot_param} is used to plot the main effect parameters fitted LSIRM model. 
#' 
#' @param data matrix; binary item response matrix to be analyzed.
#' @param lsrm_result List; The output list obtained by any lsrm function.
#' @param option character; If value is "beta", draw the boxplot for the posterior samples of beta. If value is "theta", draw the distribution of the theta estimates per total test score for the \code{data}. If value is "alpha", draw the boxplot for the posterior samples of alpha. The "alpha" is only available for 2pl LSIRM.
#' @param missing Numeric; a number to replace missing values. default value is 99.
#' 
#' @return \code{plot_param} returns the box plot of main effect parameters \eqn{\beta_i} and \eqn{\theta_j}. For item effect \eqn{\beta_i}, it shows the 95\% posterior credible intervals and for respondent effect \eqn{\theta_j}, it shows the distribution of the estimates per total sum of positive response.
#' 
#' @examples 
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm1pl(data = data)
#' plot_param(data, lsirm_result, "theta")
#' plot_param(data, lsirm_result, "beta")
#' }
#' @export
plot_param = function(data, lsrm_result, option, missing = 99){
  x <- NULL; value <- NULL
  if(option == "beta"){
    beta_dataframe <- data.frame(x = rep(1:ncol(lsrm_result$beta), each= nrow(lsrm_result$beta)),
                                 value = as.vector(lsrm_result$beta))
    #outlier 
    lower <- min(data.frame(beta_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
    upper <- max(data.frame(beta_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5, 2])
    ylim1 = c(lower, upper)
    
    if(ncol(data) > 30){
      ggplot(data=beta_dataframe, aes(x=x,y=value, group=x)) + geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = ylim1*1.05) +
        scale_x_continuous(breaks = round(seq(from = 0, to = ncol(lsrm_result$beta), length.out = 10))) +
        xlab("Item number") + ylab("Beta estimates")+
        theme(axis.text.x = element_text(face="bold",size=13),
              axis.text.y = element_text( face="bold",size=15),
              axis.title=element_text(size=15, face='bold'))
    }else{
      ggplot(data=beta_dataframe, aes(x=x,y=value, group=x))+ geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = ylim1*1.05) +
        scale_x_continuous(breaks = 0:ncol(lsrm_result$beta)) +
        xlab("Item number") + ylab("Beta estimates")+
        theme(axis.text.x = element_text(face="bold",size=13),
              axis.text.y = element_text( face="bold",size=15),
              axis.title=element_text(size=15, face='bold'))
    }
  }else if(option == "theta"){
    
    data[data==missing] = NA
    total_score <- rowSums(data, na.rm = TRUE)
    
    if(sum(data != 1 & data != 0, na.rm = TRUE) > 1){
      # Continuous Data
      binning <- cut(total_score, breaks = seq(from = min(total_score), to = max(total_score),
                                               length.out = 11), include.lowest = TRUE)
      
      theta_temp <- data.frame(x = binning, value = lsrm_result$theta_estimate)
      
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
      theta_temp <- data.frame(x = total_score, value = lsrm_result$theta_estimate)
      if(ncol(data) > 30){
        ggplot(data=theta_temp, aes(x=x,y=value, group=x))+ geom_boxplot() +
          xlab("Sum score") + ylab("Theta estimates")+
          scale_x_continuous(breaks = round(seq(from = 0, to = ncol(lsrm_result$beta), length.out = 10))) +
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold')) 
      }else{
        ggplot(data=theta_temp, aes(x=x, y=value, group = x))+ geom_boxplot() +
          xlab("Sum score") + ylab("Theta estimates")+
          scale_x_continuous(breaks = 0:ncol(lsrm_result$beta)) +
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold')) 
      }
    }
  }else if(option == "alpha"){
    if(is.null(lsrm_result$alpha) == TRUE){
      stop('The option "alpha" is only available for 2pl LSIRM.')
    }else{
      alpha_dataframe <- data.frame(x = rep(1:ncol(lsrm_result$alpha), each= nrow(lsrm_result$alpha)),
                                    value = as.vector(lsrm_result$alpha))
      
      #outlier 
      lower <- min(data.frame(alpha_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
      upper <- max(data.frame(alpha_dataframe %>% group_by(x) %>% summarise(stat = boxplot.stats(value)$stats, .groups = 'drop'))[1:ncol(data)*5, 2])
      ylim1 = c(lower, upper)
      
      if(ncol(data) > 30){
        ggplot(data=alpha_dataframe, aes(x=x,y=value, group=x)) + geom_boxplot(outlier.shape = NA) +
          coord_cartesian(ylim = ylim1*1.05) +
          scale_x_continuous(breaks = round(seq(from = 0, to = ncol(lsrm_result$alpha), length.out = 10))) +
          xlab("Item number") + ylab("Alpha estimates")+
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold'))
      }else{
        ggplot(data=alpha_dataframe, aes(x=x,y=value, group=x))+ geom_boxplot(outlier.shape = NA) +
          coord_cartesian(ylim = ylim1*1.05) +
          scale_x_continuous(breaks = 0:ncol(lsrm_result$alpha)) +
          xlab("Item number") + ylab("Alpha estimates")+
          theme(axis.text.x = element_text(face="bold",size=13),
                axis.text.y = element_text( face="bold",size=15),
                axis.title=element_text(size=15, face='bold'))
      }
    }
  }
}
