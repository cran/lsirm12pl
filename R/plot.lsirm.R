#' Plotting the interaction map or summarizing the parameter estimate of fitted LSIRM with box plot.
#'
#' @description \link{plot} is used to plot the interaction map of fitted LSIRM or summarizing the parameter estimate of fitted LSIRM with box plot.
#'
#' @param object Object of class \code{lsirm}.
#' @param option Character; If value is "interaction", draw the interaction map that represents interactions between respondents and items. If value is "beta", draw the boxplot for the posterior samples of beta. If value is "theta", draw the distribution of the theta estimates per total test score for the \code{data}. If value is "alpha", draw the boxplot for the posterior samples of alpha. The "alpha" is only available for 2PL LSIRM.
#' @param rotation Logical; If TRUE the latent positions are visualized after oblique (oblimin) rotation.
#' @param cluster Character; If value is "neyman" the cluster result are visualized by Point Process Cluster Analysis. If value is "spectral", spectral clustering method applied. Default is NA.
#' @param which.clust Character; Choose which values to clustering. "resp" is the option for respondent and "item" is the option for items. Default is "item".
#' @param interact Logical; If TRUE, draw the interaction map interactively.
#' @param chain.idx Numeric; Index of MCMC chain. Default is 1.
#' @param ... Additional arguments for the corresponding function.
#'
#' @return \code{plot} returns the interaction map or boxplot for parameter estimate.
#' @examples
#' \donttest{
#'
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5), ncol=10, nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl())
#' plot(lsirm_result)
#'
#' # use oblique rotation
#' plot(lsirm_result, rotation = TRUE)
#'
#' # interaction map interactively
#' plot(lsirm_result, interact = TRUE)
#'
#' # clustering the respondents or items
#' plot(lsirm_result, cluster = TRUE)
#' }
#' @export
plot <- function(object, ..., option = "interaction", rotation=FALSE, cluster=NA,
                 which.clust="item", interact=FALSE, chain.idx = 1){
  UseMethod("plot")
}

#' @export
plot.lsirm <- function(object, ..., option = "interaction", rotation=FALSE, cluster=NA,
                       which.clust="item", interact=FALSE, chain.idx = 1){

  group <- NULL
  type <- NULL
  value <- NULL
  if(object$chains == 1){
    x = object
  }else{
    if(chain.idx %in% 1:object$chains){
      x = object[[chain.idx]]
    }else{
      stop(sprintf("Invalid chain index: %d. The index of a chain must be lower than the total number of chains, which is %d.", chain.idx, object$chains))
    }
  }

  if(option == "interaction"){

    item_position = x$w_estimate
    resp_position = x$z_estimate
    notation = c('w','z')

    if((dim(item_position)[2]!=2)|(dim(resp_position)[2]!=2)){
      stop('\"interaction \" option is implemented for two-dimensional latent space.')
    }

    if(is.na(cluster)){
      axis1 <- NULL; axis2 <- NULL

      if(rotation){
        rot <- oblimin(item_position)
        item_position = rot$loadings
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


      # plotly
      if(interact){
        df1$type = rep("item", nrow(df1))
        df2$type = rep("res", nrow(df2))
        df.interact=rbind(df1,df2)
        colnames(df.interact)=c('axis1','axis2','source','type')
        pos.label = rep(c("item","res"),
                        times=c(nrow(df1),nrow(df2)))
        label.n = paste0(pos.label,
                         c(1:nrow(df1),1:nrow(df2)))
        item = c(1:nrow(df1), rep("",nrow(df2)))
        interact = ggplot(data = df.interact,
                          aes(x=axis1, y=axis2, color=type,
                              lb=label.n,
                              text=paste0("x.axis: ",round(axis1,3), "<br>",
                                          "y.axis: ",round(axis2,3), "<br>")))+
          geom_point(size = ifelse(df.interact$type=="item",
                                   1e-06,1)) +
          geom_text(data = df.interact, aes(x=axis1, y=axis2, label=item))+
          scale_color_manual(values=c("res"="black","item"="red"))+
          xlim(axis_range) +
          ylim(axis_range) +
          coord_cartesian(expand = FALSE) + theme_bw() +
          theme(plot.margin = unit(c(1,1,1,1), "cm"),
                axis.text=element_text(size=16),
                axis.title=element_text(size=14,face="bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.title=element_blank(),
                # legend.position.inside = c(0.9,0.9),
                legend.text = element_text(size=10),
                plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
          ggtitle("Interaction Map")
        ggplotly(interact, tooltip = c("text","lb"))
      }else{
        ggplot() +
          geom_point(data = df2, aes(x = axis1, y = axis2), size = 0.7) +
          geom_text(data = df1, aes(x = axis1, y = axis2, label=1:nrow(df1)),
                    color = "red", size = 4, fontface = "bold") +
          xlim(axis_range)+ylim(axis_range) + coord_cartesian(expand = FALSE) + theme_bw() +
          theme(plot.margin = unit(c(1,1,1,1), "cm"),
                axis.text=element_text(size=16),
                axis.title=element_text(size=14,face="bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.title=element_blank(),
                # legend.position.inside = c(0.9,0.9),
                legend.text = element_text(size=16),
                plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
          ggtitle("Interaction Map")
      }
    }else{


      if(cluster == "neyman"){
        if(which.clust == "item"){
          w.samps1 <- item_position
          z.samps1 <- resp_position
        }else if(which.clust == "resp"){
          w.samps1 <- resp_position
          z.samps1 <- item_position
        }else{
          stop("Not supported")
        }

        df = rbind(item_position, resp_position)
        max_coord = apply(df,2,max,na.rm=T)
        min_coord = apply(df,2,min,na.rm=T)

        W <- owin(xrange=c(0,1), yrange=c(0,1)) # spatstat package ()

        # Normalizing the w
        # x <- (w.samps1[,1] - min(w.samps1[,1])) / (max(w.samps1[,1]) - min(w.samps1[,1]))
        # y <- (w.samps1[,2] - min(w.samps1[,2])) / (max(w.samps1[,2]) - min(w.samps1[,2]))
        x = (w.samps1[,1] - min_coord[1])/(max_coord[1]-min_coord[1])
        y = (w.samps1[,2] - min_coord[2])/(max_coord[2]-min_coord[2])
        w.post1 <- data.frame(cbind(x, y))
        colnames(w.post1) <- c('x', 'y')

        # Normalizing the z
        # zx <- (z.samps1[,1] - min(z.samps1[,1])) / (max(z.samps1[,1]) - min(z.samps1[,1]))
        # zy <- (z.samps1[,2] - min(z.samps1[,2])) / (max(z.samps1[,2]) - min(z.samps1[,2]))
        zx = (z.samps1[,1] - min_coord[1])/(max_coord[1]-min_coord[1])
        zy = (z.samps1[,2] - min_coord[2])/(max_coord[2]-min_coord[2])
        z.post1 <- data.frame(cbind(zx, zy))
        colnames(z.post1) <- c('x', 'y')

        xlim <- W$xrange; ylim <- W$yrange
        AreaW <- 1 # |S| area of the interaction map domain

        U <- w.post1
        x <- U[,1]
        y <- U[,2]
        X <- t(rbind(x, y))

        # alpha: an expected number of items for each group center ci
        # omega: controls the range of item groups in the interaction map
        salpha <- 0.1; somega <- 0.015
        sd_alpha <- 0.1; sd_omega <- 0.015
        Niter = 100

        parent <- list()
        parentnum <- c()
        accept <- c()
        logllh <- c()
        alphalist <- list()
        kappalist <- list()
        omegalist <- list()

        pb <- txtProgressBar(title = "progress bar", min = 0, max = Niter,
                             style = 3, width = 50)


        # To calculate the BIC repeat 100 times
        for (i in 1:Niter){
          setTxtProgressBar(pb, i, label = paste( round(i/Niter * 100, 0), "% done"))

          # Update the alpha, omega, CC (Thomas process fitting procedures)
          Thomas<-MCMCestThomas(X, xlim, ylim, NStep=25000, DiscardStep=5000, Jump=5)

          # Summarize the result of Thomas process fitting process
          CC <- Thomas$CC
          omega <- Thomas$Omegahat
          alpha <- Thomas$Alphahat
          integral <- Kumulavsech(CC, omega, xlim, ylim)
          logP <- logpXCbeta(X, CC, alpha, omega, AreaW, integral)

          # Save the result for each repeat
          parentnum <- c(parentnum, dim(Thomas$CC)[1])
          parent[[i]] <- CC
          accept <- c(accept, Thomas$Accept)
          logllh <- c(logllh, logP)
          alphalist[[i]] <- Thomas$Postalpha
          kappalist[[i]] <- Thomas$Postkappa
          omegalist[[i]] <- Thomas$Postomega
        }

        # BIC
        bic.total <- bic(X, parentnum, logllh)
        df <- data.frame(parentnum, bic.total)

        ## BIC of the most frequent parent number
        ind <- as.numeric(names(which.max(table(parentnum))))

        ## The index number which parent number is "ind (most frequent parent number)"
        ind5 <- c()
        for(i in 1:Niter){
          if(parentnum[i]==ind){
            ind5 <- c(ind5, i)
          }
        }

        bic5 <- bic.total[ind5]

        # The index number having smallest BIC among number of parent is "ind"
        ind2 <- as.numeric(rownames(df[((df["parentnum"]==ind) & df["bic.total"]==bic5[which.min(bic5)] ),]))

        # location of "in2" (which is best interation whith having smallest BIC)
        par5 <- parent[[ind2]]
        par5 <- data.frame(par5)
        colnames(par5) <- c('x', 'y')

        # clustering the latent position of item (using distance)
        g_item <- clustering(w.post1, par5)
        colnames(g_item) <- c("distance", "group")
        # raanking based on distance from center
        # (The closer the center is, the greater the alpha is.)
        g_alpha <- c()

        for(l in 1:ind){
          temp <- g_item[which(g_item$group==l),]
          temp <- ranking(temp)
          g_alpha <- rbind(g_alpha, temp)
        }

        # ordering by order of item
        g_alpha <- g_alpha[order(as.numeric(rownames(g_alpha))),]
        g_fin <- cbind(w.post1[,1:2], g_alpha)

        # parent density using the 1000 estimated parent location
        cc <- parent[[1]]
        for(i in 2:Niter){
          cc<-rbind(cc, parent[[i]])
        }
        cc <- data.frame(cc)
        colnames(cc) <- c('x','y')

        if(ind > length(LETTERS)){
          addlabel = LETTERS
          for(U in 1:length(LETTERS)){
            for(l in 1:length(letters)){
              addlabel = c(addlabel, paste0(LETTERS[U],letters[l]))
            }
          }
          alphabet = addlabel[1:ind]
        } else{
          alphabet = LETTERS[1:ind]
        }

        # Set the color and cluster name
        ggcolor = rainbow(ind, s=.6, v=.9)[sample(1:ind, ind)]
        par5["cluster"] <- alphabet
        par5["color"] <- ggcolor


        # plotly
        if(interact){
          temp <- (ggplot(w.post1, aes(x, y)) +
                     geom_point(data=z.post1, aes(x, y), col="grey", cex=1.0) +
                     stat_density_2d(data=cc, aes(x, y), color="gray80") + #density
                     geom_text(data=g_fin, aes(x, y), label=rownames(g_fin),
                               color=ggcolor[g_fin$group], cex=4, fontface="bold") + # number of item
                     geom_text(data=par5, aes(x, y), label=alphabet[1:ind], col="gray30", cex=5.5, fontface="bold") + #alphabet
                     theme_bw() +
                     theme(plot.margin = unit(c(1,1,1,1), "cm"),
                           axis.text=element_text(size=16),
                           axis.title=element_text(size=14,face="bold"),
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank(),
                           legend.title=element_blank(),
                           legend.text = element_text(size=16),
                           plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
                     ggtitle("Interaction Map"))
          if(which.clust == "item"){
            int.plot = ggplotly(temp) %>%
              add_markers(x = z.post1$x, y = z.post1$y,
                          type = 'scatter',
                          mode = 'markers',
                          text = paste("respondent", 1:nrow(z.post1), sep = ""),
                          name = "Respondent",
                          marker = list(color = "lightgrey")) %>%
              add_text(x = g_fin$x, y = g_fin$y, type = 'scatter',
                       name = "Item",
                       text = 1:nrow(w.post1),
                       textfont = list(family="Arial Black",size=16, weight="bold", color = ggcolor[g_fin$group])) %>%
              add_text(x = par5$x, y = par5$y, type = 'scatter',
                       name = "Cluster",
                       text = alphabet[1:ind],
                       textfont = list(size=19, color = "black"))
          }else if(which.clust == "resp"){
            int.plot = ggplotly(temp) %>%
              add_markers(x = z.post1$x, y = z.post1$y,
                          type = 'scatter',
                          mode = 'markers',
                          text = paste("item", 1:nrow(z.post1), sep = ""),
                          name = "Item",
                          marker = list(color = "lightgrey")) %>%
              add_text(x = g_fin$x, y = g_fin$y, type = 'scatter',
                       name = "Respondent",
                       text = 1:nrow(w.post1),
                       textfont = list(family="Arial Black",size=16, weight="bold", color = ggcolor[g_fin$group])) %>%
              add_text(x = par5$x, y = par5$y, type = 'scatter',
                       name = "Cluster",
                       text = alphabet[1:ind],
                       textfont = list(size=19, color = "black"))
          }else{
            stop("Not supported")
          }

          print(int.plot)

        }else{
          # Draw plot using ggplot
          print(ggplot(w.post1, aes(x, y)) +
                  geom_point(data=z.post1, aes(x, y), col="grey", cex=1.0) +
                  stat_density_2d(data=cc, aes(x, y), color="gray80") + #density
                  geom_text(data=g_fin, aes(x, y), label=rownames(g_fin), color=ggcolor[g_fin$group], cex=4, fontface="bold") + # number of item
                  geom_text(data=par5, aes(x, y), label=alphabet[1:ind], col="gray30", cex=5.5, fontface="bold") + #alphabet
                  theme_bw() +
                  theme(plot.margin = unit(c(1,1,1,1), "cm"),
                        axis.text=element_text(size=16),
                        axis.title=element_text(size=14,face="bold"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.title=element_blank(),
                        legend.text = element_text(size=16),
                        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
                  ggtitle("Interaction Map"))
        }

        # print the cluster information
        g_alpha <- cbind(g_alpha, item = rownames(g_alpha))
        clust <- data.frame(
          g_alpha %>%
            group_by(group) %>%
            reframe(items = paste(item, collapse = ", ")))



        c_res = capture.output(print(data.frame(group = alphabet[clust$group], item = clust[,2]),
                                     row.names=F))
        cat("\n\nClustering result (Neyman-Scott process): \n",paste(c_res,"\n",sep=" "))

      }else if(cluster == "spectral"){

        # Select the number of clustering using Average Silhouette Width (ASW) which is a popular cluster validation index to estimate the number of clusters.

        scale_func <- function(x) x %>% mutate_if(is.numeric, function(y) as.vector(scale(y)))
        if(which.clust == "item"){
          pos_scaled <- item_position %>% scale()
          df1=as.data.frame(item_position)
          df2=as.data.frame(resp_position)
        }else if(which.clust == "resp"){
          pos_scaled <- resp_position %>% scale()
          df1=as.data.frame(resp_position)
          df2=as.data.frame(item_position)
        }else{
          stop("Not supported")
        }

        num <- 2:(nrow(pos_scaled)-1)
        ASW <- sapply(num, FUN=function(k) {
          fpc::cluster.stats(dist(pos_scaled),
                             kmeans(pos_scaled, centers=k,
                                    nstart = 5)$cluster)$avg.silwidth
        })

        best_k <- num[which.max(ASW)]
        # asw = pbapply::pbsapply(num, FUN=function(k) {
        #   fpc::cluster.stats(dist(pos_scaled),
        #                      kmeans(pos_scaled, centers=k,
        #                             nstart = 5)$cluster)$avg.silwidth
        # })
        # Spectral Clustering
        spectral_result <- specc(as.matrix(pos_scaled), centers = best_k) #kernlab package
        spectral_result <- data.frame(cbind(group = as.numeric(spectral_result), item = 1:nrow(pos_scaled)))
        clust <- data.frame(
          spectral_result %>%
            group_by(group) %>%
            reframe(items = paste(item, collapse = ", "))
          # summarise(items = paste(item, collapse = ", "))
        )

        df1[,3]=notation[1]
        df2[,3]=notation[2]

        colnames(df1)=c('x','y','source')
        colnames(df2)=c('x','y','source')


        df=rbind(df2,df1)
        colnames(df)=c('axis1','axis2','source')

        max_coordinate = sapply(df[,c(1,2)], max, na.rm = TRUE)
        min_coordinate = sapply(df[,c(1,2)], min, na.rm = TRUE)
        axis_value = max(abs(c(max_coordinate,min_coordinate)))
        axis_range = c(-axis_value,axis_value)*1.1

        g_fin <- cbind(df1, spectral_result)
        ind <- max(spectral_result$group)

        if(ind > length(LETTERS)){
          addlabel = LETTERS
          for(U in 1:length(LETTERS)){
            for(l in 1:length(letters)){
              addlabel = c(addlabel, paste0(LETTERS[U],letters[l]))
            }
          }
          alphabet = addlabel[1:ind]
        } else{
          alphabet = LETTERS[1:ind]
        }

        # Set the color and cluster name
        ggcolor = rainbow(ind, s=.6, v=.9)[sample(1:ind, ind)]

        # plotly
        if(interact){
          temp <- (ggplot(data=df, aes(x, y)) +
                     geom_point(data=df2, aes(x, y), col="grey", cex=1.0) +
                     geom_text(data=g_fin, aes(x, y), label=g_fin$item,
                               color=ggcolor[g_fin$group], cex=4, fontface="bold") + # number of item
                     theme_bw() +
                     theme(plot.margin = unit(c(1,1,1,1), "cm"),
                           axis.text=element_text(size=16),
                           axis.title=element_text(size=14,face="bold"),
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank(),
                           legend.title=element_blank(),
                           legend.text = element_text(size=16),
                           plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
                     ggtitle("Interaction Map"))
          if(which.clust == "item"){
            int.plot = ggplotly(temp) %>%
              add_markers(x = df2$x, y = df2$y,
                          type = 'scatter',
                          mode = 'markers',
                          text = paste("respondent", 1:nrow(df2), sep = ""),
                          name = "Respondent",
                          marker = list(color = "lightgrey")) %>%
              add_text(x = g_fin$x, y = g_fin$y, type = 'scatter',
                       name = "Item",
                       text = 1:nrow(g_fin),
                       textfont = list(family="Arial Black",size=16, weight="bold", color = ggcolor[g_fin$group]))
          }else if(which.clust == "resp"){
            int.plot = ggplotly(temp) %>%
              add_markers(x = df2$x, y = df2$y,
                          type = 'scatter',
                          mode = 'markers',
                          text = paste("item", 1:nrow(df2), sep = ""),
                          name = "Item",
                          marker = list(color = "lightgrey")) %>%
              add_text(x = g_fin$x, y = g_fin$y, type = 'scatter',
                       name = "Respondent",
                       text = 1:nrow(g_fin),
                       textfont = list(family="Arial Black",size=16, weight="bold", color = ggcolor[g_fin$group]))
          }else{
            stop("Not supported")
          }
          print(int.plot)
        }else{
          temp <- (ggplot(data=df, aes(x, y)) +
                     geom_point(data=df2, aes(x, y), col="grey", cex=1.0) +
                     geom_text(data=g_fin, aes(x, y), label=g_fin$item,
                               color=ggcolor[g_fin$group], cex=4, fontface="bold") + # number of item
                     theme_bw() +
                     theme(plot.margin = unit(c(1,1,1,1), "cm"),
                           axis.text=element_text(size=16),
                           axis.title=element_text(size=14,face="bold"),
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank(),
                           legend.title=element_blank(),
                           legend.text = element_text(size=16),
                           plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
                     ggtitle("Interaction Map"))
          print(temp)
        }
        # print the cluster information

        c_res = capture.output(print(data.frame(group = alphabet[clust$group], item = clust[,2]),
                                     row.names=F))
        cat("\n\nClustering result (Spectral Clustering): \n",paste(c_res,"\n",sep=" "))


      }
    }
  }else if(option == "beta"){

    if(option == "beta" & !is.na(cluster)){
      stop('\"cluster\" option is implemented for interaction map.')
    }

    beta_dataframe <- data.frame(x = rep(1:ncol(x$beta), each= nrow(x$beta)),
                                 value = as.vector(x$beta))
    #outlier
    lower <- min(data.frame(beta_dataframe %>% group_by(x) %>%
                              reframe(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
    # summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
    upper <- max(data.frame(beta_dataframe %>% group_by(x) %>%
                              reframe(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5, 2])
    # summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5, 2])
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

    if(option == "theta" & !is.na(cluster)){
      stop('\"cluster\" option is implemented for interaction map.')
    }


    if(is.null(x$missing.val) == TRUE){
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
    }else if(option == "alpha" & !is.na(cluster)){
      stop('\"cluster\" option is implemented for interaction map.')
    }else{
      alpha_dataframe <- data.frame(x = rep(1:ncol(x$alpha), each= nrow(x$alpha)),
                                    value = as.vector(x$alpha))

      #outlier
      lower <- min(data.frame(alpha_dataframe %>% group_by(x) %>%
                                reframe(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
      # summarise(stat = boxplot.stats(value)$stats, .groups = "drop"))[1:ncol(data)*5 - 4, 2])
      upper <- max(data.frame(alpha_dataframe %>% group_by(x) %>%
                                reframe(stat = boxplot.stats(value)$stats, .groups = 'drop'))[1:ncol(data)*5, 2])
      # summarise(stat = boxplot.stats(value)$stats, .groups = 'drop'))[1:ncol(data)*5, 2])
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
