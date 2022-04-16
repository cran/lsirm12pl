#' Plotting the latent space of fitted LSIRM model
#' 
#' @description \link{plot_latent} is used to plot the latent space of fitted LSIRM model.
#' 
#' @param lsrm_result List; The output list obtained by any lsrm function.
#' @param rotation Logical; If TRUE the latent positions are visualized after oblique (oblimin) rotation.
#' 
#' @return \code{plot_latent} returns the plot of latent space visualize an interaction map that represents interactions between respondents and items. 
#' 
#' @examples 
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm1pl(data = data)
#' plot_latent(lsirm_result)
#' # use oblique rotation
#' plot_latent(lsirm_result, rotation = TRUE) 
#' }
#' @export
plot_latent = function(lsrm_result, rotation=FALSE) {
  axis1 <- NULL; axis2 <- NULL
  # globalVariables(c("axis1", "axis2"))
  if(sum(lsrm_result$z_estimate) == 0){
    item_position = lsrm_result$ls_mean_item
    resp_position = lsrm_result$ls_mean_respondent
    notation = c('v','u')
  }else{
    item_position = lsrm_result$w_estimate
    resp_position = lsrm_result$z_estimate
    notation = c('w','z')
  }
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
  
  # ggplot(df,aes(axis1,axis2,col=source)) + 
  #   xlim(axis_range)+ylim(axis_range) + coord_cartesian(expand = FALSE) + theme_bw() + 
  #   geom_point(size = 2,alpha=1,shape=16) + 
  #   theme(plot.margin = unit(c(1,1,1,1), "cm"),
  #         axis.text=element_text(size=16),
  #         axis.title=element_text(size=14,face="bold"),
  #         axis.title.x=element_blank(),
  #         axis.title.y=element_blank(),
  #         legend.title=element_blank(),
  #         legend.position = c(0.9,0.9),
  #         legend.text = element_text(size=16)
  #   )
  
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
          legend.text = element_text(size=16)
    )
  
}

