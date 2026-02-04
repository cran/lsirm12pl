#' Diagnostic the result of LSIRM.
#'
#' @description \code{diagnostic} checks the convergence of MCMC for LSIRM parameters using various diagnostic tools, such as trace plots, posterior density distributions, autocorrelation functions (ACF), and Gelman-Rubin-Brooks plots.
#'
#' @param object Object of class \code{lsirm}.
#' @param draw.item List; Each key in the list corresponds to a specific parameters such as "beta", "theta", "gamma", "alpha", "sigma", "theta_sd", "z", "w", and "zw.dist". The values of the list indicate the indices of these parameters. For the key "zw.dist", the value is a matrix with two columns: the first column represents the indices of respondents, and the second column represents the indices of items. For the keys "z" and "w", the value is a matrix with two columns: the first column represents the indices of respondents/items, and the second column represents the dimension indices.
#' @param gelman.diag Logical; If TRUE, the Gelman-Rubin convergence diagnostic will be printed. Default is FALSE.
#'
#' @return \code{diagnostic} returns plots for checking MCMC convergence for selected parameters.
#'
#' @examples
#' \donttest{
#' # Generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5), ncol=10, nrow=50)
#'
#' # For 1PL LSIRM
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = FALSE, fixed_gamma = FALSE))
#' diagnostic(lsirm_result)
#'
#' # For 2PL LSIRM
#' lsirm_result <- lsirm(data ~ lsirm2pl(spikenslab = FALSE, fixed_gamma = FALSE))
#' diagnostic(lsirm_result)
#'
#' # For specific latent positions
#' diagnostic(lsirm_result, draw.item = list(z = matrix(c(1,1, 2,1, 1,2), ncol=2, byrow=TRUE),
#'                                        w = matrix(c(1,1, 2,1), ncol=2, byrow=TRUE)))
#'
#' }
#' @export diagnostic
diagnostic <- function(object,
                       draw.item = list(beta=c(1),
                                        theta=c(1)),
                       gelman.diag = FALSE){
  UseMethod("diagnostic")
}

#' @export
diagnostic.lsirm <- function(object,
                             draw.item = list(beta=c(1),
                                              theta=c(1)),
                             gelman.diag = FALSE)
{

  ACF <- Chain <- Iteration <- Lag <- PSRF <- Type <- iteration <- var1 <- NULL
  orders = data.frame(idx = c(1:9),
                      param = c("beta", "theta", "gamma", "alpha", "sigma", "theta_sd", "z", "w", "zw.dist"))
  which.draw = names(draw.item)
  porder = orders[orders$param %in% which.draw,]

  if(object$chains > 1){
    multi_chain = TRUE
    nmcmc <- dim(object[[1]]$w)[1]
    which.draw.temp = which.draw
  }else{
    multi_chain = FALSE
    nmcmc <- dim(object$w)[1]
    which.draw.temp = which.draw
  }
  if(multi_chain){

    if((is.null(object[[1]]$beta)) & ("beta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$gamma)) & ("gamma" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$theta)) & ("theta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$z) | is.null(object[[1]]$w)) & ("zw.dist" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$alpha)) & ("alpha" %in% which.draw)) stop("MCMC sample was not stored. The parameter alpha is only stored in method lsirm2pl")

  }else{
    if((is.null(object$beta)) & ("beta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$gamma)) & ("gamma" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$theta)) & ("theta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$z) | is.null(object$w)) & ("zw.dist" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$alpha)) & ("alpha" %in% which.draw)) stop("MCMC sample was not stored. The parameter alpha is only stored in method lsirm2pl")

  }

  if(is.null(draw.item$beta)&("beta" %in% which.draw)) draw.item$beta = "first"
  if(is.null(draw.item$theta)&("theta" %in% which.draw)) draw.item$theta = "first"
  if(is.null(draw.item$alpha)&("alpha" %in% which.draw)) draw.item$alpha = 2
  if(is.null(draw.item$z)&("z" %in% which.draw)) draw.item$z = matrix(c(1, 1), ncol = 2)
  if(is.null(draw.item$w)&("w" %in% which.draw)) draw.item$w = matrix(c(1, 1), ncol = 2)
  if(is.null(draw.item$zw.dist)&("zw.dist" %in% which.draw)) draw.item$zw.dist = matrix(c(1, 1), ncol = 2)
  chain_list_all <- vector("list", length = object$chains)

  .gelman_plot_safe <- function(mcmclist, param){
    mcmclist_param <- coda::as.mcmc.list(lapply(mcmclist, function(ch){
      mat <- as.matrix(ch)
      coda::mcmc(mat[, param, drop = FALSE], thin = 1)
    }))

    # Open a null device to capture plot output without displaying
    current_dev <- grDevices::dev.cur()
    grDevices::pdf(file = NULL)
    
    # Try to get gelman.plot data
    out <- tryCatch({
      coda::gelman.plot(mcmclist_param, auto.layout = FALSE)
    }, error = function(e){
      grDevices::dev.off()
      if(current_dev > 1) grDevices::dev.set(current_dev)
      message("Skipping Gelman-Rubin plot: ", conditionMessage(e))
      array_stub <- array(
        NA_real_,
        dim = c(1, 1, 2),
        dimnames = list(NULL, NULL, c("median", "97.5%"))
      )
      return(list(shrink = array_stub, last.iter = NA_integer_))
    })
    
    # Close the null device and restore previous device
    grDevices::dev.off()
    if(current_dev > 1) grDevices::dev.set(current_dev)

    out
  }

  readline <- function(prompt = ""){
    if(interactive()) base::readline(prompt = prompt) else ""
  }

  #### Beta  ------------
  if("beta" %in% which.draw){

    if(multi_chain){
      is_grm <- !is.null(object[[1]]$method) && object[[1]]$method == "lsirmgrm"
      
      # Handle GRM beta as list structure for each chain
      for(i in 1:object$chains){
        if(is_grm && is.list(object[[i]]$beta) && !is.matrix(object[[i]]$beta)){
          # Beta is stored as list with threshold components (th1, th2, etc.)
          beta_list <- object[[i]]$beta
          nitem <- ncol(object[[i]]$data)
          nthreshold <- length(beta_list)
          item_names <- if(!is.null(colnames(object[[i]]$data))) colnames(object[[i]]$data) else paste0("i", 1:nitem)
          
          # Combine all thresholds by column
          beta_matrix <- do.call(cbind, beta_list)
          
          # Create column names: item:threshold format
          threshold_names <- names(beta_list)
          colnames_beta <- character(nitem * nthreshold)
          for(item in 1:nitem){
            for(th in 1:nthreshold){
              colnames_beta[(item-1)*nthreshold + th] <- paste0(item_names[item], ":", threshold_names[th])
            }
          }
          colnames(beta_matrix) <- colnames_beta
          object[[i]]$beta <- beta_matrix
        }
      }
      
      if(is.null(colnames(object[[1]]$beta))){
        colnames(object[[1]]$beta) = 1:ncol(object[[1]]$beta)
      }else{
        if(!is_grm && !is.null(colnames(object[[1]]$data))){
          colnames(object[[1]]$beta) <- colnames(object[[1]]$data)
        }
      }



      if(is_grm && is.matrix(draw.item$beta) && ncol(draw.item$beta) == 2){
        Kminus1 <- object[[1]]$ncat - 1
        item_idx <- draw.item$beta[,1]
        th_idx <- draw.item$beta[,2]
        if(any(item_idx < 1 | item_idx > ncol(object[[1]]$data))) stop("Invalid item index in draw.item$beta")
        if(any(th_idx < 1 | th_idx > Kminus1)) stop("Invalid threshold index in draw.item$beta")
        draw.item.num <- (item_idx - 1) * Kminus1 + th_idx
        draw.item.temp <- colnames(object[[1]]$beta)[draw.item.num]
      }else if((length(draw.item$beta) == 1) & (draw.item$beta[1] == "first")){
        draw.item.temp = colnames(object[[1]]$beta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$beta)){
        draw.item.temp = colnames(object[[1]]$beta)[draw.item$beta]
        draw.item.num <- draw.item$beta
      }else{
        draw.item.temp = draw.item$beta
        draw.item.num = which(colnames(object[[1]]$beta) %in% draw.item.temp)
      }


      pnames <- c(paste('beta [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$beta[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{
      is_grm <- !is.null(object$method) && object$method == "lsirmgrm"
      
      # Handle GRM beta as list structure
      if(is_grm && is.list(object$beta) && !is.matrix(object$beta)){
        # Beta is stored as list with threshold components (th1, th2, etc.)
        # Each component is a matrix: (n_samples x n_items)
        # Convert to single matrix: (n_samples x (n_items * n_thresholds))
        beta_list <- object$beta
        nitem <- ncol(object$data)
        nthreshold <- length(beta_list)
        item_names <- if(!is.null(colnames(object$data))) colnames(object$data) else paste0("i", 1:nitem)
        
        # Combine all thresholds by column
        beta_matrix <- do.call(cbind, beta_list)
        
        # Create column names: item:threshold format
        threshold_names <- names(beta_list)
        colnames_beta <- character(nitem * nthreshold)
        for(i in 1:nitem){
          for(j in 1:nthreshold){
            colnames_beta[(i-1)*nthreshold + j] <- paste0(item_names[i], ":", threshold_names[j])
          }
        }
        colnames(beta_matrix) <- colnames_beta
        object$beta <- beta_matrix
      }
      
      if(is.null(colnames(object$beta))){
        colnames(object$beta) = 1:ncol(object$beta)
      }else{
        if(!is_grm && !is.null(colnames(object$data))){
          colnames(object$beta) <- colnames(object$data)
        }
      }



      if(is_grm && is.matrix(draw.item$beta) && ncol(draw.item$beta) == 2){
        Kminus1 <- object$ncat - 1
        item_idx <- draw.item$beta[,1]
        th_idx <- draw.item$beta[,2]
        if(any(item_idx < 1 | item_idx > ncol(object$data))) stop("Invalid item index in draw.item$beta")
        if(any(th_idx < 1 | th_idx > Kminus1)) stop("Invalid threshold index in draw.item$beta")
        draw.item.num <- (item_idx - 1) * Kminus1 + th_idx
        draw.item.temp <- colnames(object$beta)[draw.item.num]
      }else if((length(draw.item$beta) == 1) & (draw.item$beta[1] == "first")){
        draw.item.temp = colnames(object$beta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$beta)){
        draw.item.temp = colnames(object$beta)[draw.item$beta]
        draw.item.num <- draw.item$beta
      }else{
        draw.item.temp = draw.item$beta
        draw.item.num = which(colnames(object$beta) %in% draw.item.temp)
      }


      pnames <- c(paste('beta [', draw.item.temp, ']', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$beta[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))
      combined_data$var1 <- combined_data$value


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

          median_values <- gelman_plot$shrink[,, "median"]
          upper_values <- gelman_plot$shrink[,, "97.5%"]
          iterations <- gelman_plot$last.iter

          # Create a data frame
          psrf_data <- data.frame(
            Iteration = iterations,
            MedianPSRF = as.vector(median_values),
            UpperPSRF = as.vector(upper_values)
          )

          # Reshape for plotting
          psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
            mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

          psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
            geom_line(aes(linetype = Type, color = Type)) +
            scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
            scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
            labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
            theme_minimal() +
            theme(legend.title = element_blank(),
                  legend.position.inside = c(1, 1),
                  legend.justification = c("right", "top")

            )
        }else{
          psrf_plot <- grid::nullGrob()
        }


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "beta",]$idx == max(porder$idx)){
          if(((length(draw.item.num) > 1)&(i < length(draw.item.num))) ){
            if(interactive()) readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })



        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "beta",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("beta"))
    # coda::gelman.diag(mcmclist)
  } ## Beta

  ## Theta ----------
  if("theta" %in% which.draw){

    if(multi_chain){

      if(is.null(rownames(object[[1]]$data))){
        colnames(object[[1]]$theta) = 1:ncol(object[[1]]$theta)
      }else{colnames(object[[1]]$theta) <- rownames(object[[1]]$data)}

      if((length(draw.item$theta) == 1) & (draw.item$theta[1] == "first")){
        draw.item.temp = colnames(object[[1]]$theta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$theta)){
        draw.item.temp = colnames(object[[1]]$theta)[draw.item$theta]
        draw.item.num <- draw.item$theta
      }else{
        draw.item.temp = draw.item$theta
        draw.item.num = which(colnames(object[[1]]$theta) %in% draw.item.temp)
      }

      pnames <- c(paste('theta [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$theta[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{
      if(is.null(rownames(object$data))){
        colnames(object$theta) = 1:ncol(object$theta)
      }else{colnames(object$theta) <- rownames(object$data)}

      if((length(draw.item$theta) == 1) & (draw.item$theta[1] == "first")){
        draw.item.temp = colnames(object$theta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$theta)){
        draw.item.temp = colnames(object$theta)[draw.item$theta]
        draw.item.num <- draw.item$theta
      }else{
        draw.item.temp = draw.item$theta
        draw.item.num = which(colnames(object$theta) %in% draw.item.temp)
      }

      pnames <- c(paste('theta [', draw.item.temp, ']', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$theta[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))

    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))
      combined_data$var1 <- combined_data$value


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        gelman_plot <- .gelman_plot_safe(mcmclist, param)



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position.inside = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("theta"))
  } ## Theta

  ## Gamma =====
  if("gamma" %in% which.draw){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('gamma', sep = ''))

      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$gamma[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('gamma', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$gamma[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))
      combined_data$var1 <- combined_data$value


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

          median_values <- gelman_plot$shrink[,, "median"]
          upper_values <- gelman_plot$shrink[,, "97.5%"]
          iterations <- gelman_plot$last.iter

          # Create a data frame
          psrf_data <- data.frame(
            Iteration = iterations,
            MedianPSRF = as.vector(median_values),
            UpperPSRF = as.vector(upper_values)
          )

          # Reshape for plotting
          psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
            mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

          psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
            geom_line(aes(linetype = Type, color = Type)) +
            scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
            scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
            labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
            theme_minimal() +
            theme(legend.title = element_blank(),
                  legend.position.inside = c(1, 1),
                  legend.justification = c("right", "top")
            )
        }else{
          psrf_plot <- grid::nullGrob()
        }


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "gamma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "gamma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("gamma"))
  } ## Gamma


  ## alpha -----
  if("alpha" %in% which.draw){

    if(multi_chain){
      if(is.null(colnames(object[[1]]$data))){
        colnames(object[[1]]$alpha) = 1:ncol(object[[1]]$alpha)
      }else{colnames(object[[1]]$alpha) <- colnames(object[[1]]$data)}

      if((length(draw.item$alpha) == 1) & (draw.item$alpha[1] == "first")){
        draw.item.temp = colnames(object[[1]]$alpha)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$alpha)){
        draw.item.temp = colnames(object[[1]]$alpha)[draw.item$alpha]
        draw.item.num <- draw.item$alpha
      }else{
        draw.item.temp = draw.item$alpha
        draw.item.num = which(colnames(object[[1]]$alpha) %in% draw.item.temp)
      }

      pnames <- c(paste('alpha [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$alpha[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{
      if(is.null(colnames(object$data))){
        colnames(object$alpha) = 1:ncol(object$alpha)
      }else{colnames(object$alpha) <- colnames(object$data)}

      if((length(draw.item$alpha) == 1) & (draw.item$alpha[1] == "first")){
        draw.item.temp = colnames(object$alpha)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$alpha)){
        draw.item.temp = colnames(object$alpha)[draw.item$alpha]
        draw.item.num <- draw.item$alpha
      }else{
        draw.item.temp = draw.item$alpha
        draw.item.num = which(colnames(object$alpha) %in% draw.item.temp)
      }

      pnames <- c(paste('alpha [', draw.item.temp, ']', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$alpha[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))
      combined_data$var1 <- combined_data$value
      combined_data$var1 <- combined_data$value


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))),
                     linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

          median_values <- gelman_plot$shrink[,, "median"]
          upper_values <- gelman_plot$shrink[,, "97.5%"]
          iterations <- gelman_plot$last.iter

          # Create a data frame
          psrf_data <- data.frame(
            Iteration = iterations,
            MedianPSRF = as.vector(median_values),
            UpperPSRF = as.vector(upper_values)
          )

          # Reshape for plotting
          psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
            mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

          psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
            geom_line(aes(linetype = Type, color = Type)) +
            scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
            scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
            labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
            theme_minimal() +
            theme(legend.title = element_blank(),
                  legend.position.inside = c(1, 1),
                  legend.justification = c("right", "top")

            )
        }else{
          psrf_plot <- grid::nullGrob()
        }


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "alpha",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))
        if(porder[porder$param == "alpha",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }

  } ## Alpha


  ## sigma --------
  if("sigma" %in% which.draw){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('sigma', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        sigma_draw <- object[[i]]$sigma
        sigma_vec <- if(is.null(dim(sigma_draw))) as.numeric(sigma_draw) else as.numeric(sigma_draw[, draw.item.num])
        chain_list[[i]] <- matrix(sigma_vec, ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('sigma', sep = ''))
      chains <- 1

      chain_list <- list()
      sigma_draw <- object$sigma
      sigma_vec <- if(is.null(dim(sigma_draw))) as.numeric(sigma_draw) else as.numeric(sigma_draw[, draw.item.num])
      chain_list[[1]] <- matrix(sigma_vec, ncol = length(pnames),
                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))
      combined_data$var1 <- combined_data$value


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))),
                     linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

          median_values <- gelman_plot$shrink[,, "median"]
          upper_values <- gelman_plot$shrink[,, "97.5%"]
          iterations <- gelman_plot$last.iter

          # Create a data frame
          psrf_data <- data.frame(
            Iteration = iterations,
            MedianPSRF = as.vector(median_values),
            UpperPSRF = as.vector(upper_values)
          )

          # Reshape for plotting
          psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
            mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

          psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
            geom_line(aes(linetype = Type, color = Type)) +
            scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
            scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
            labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
            theme_minimal() +
            theme(legend.title = element_blank(),
                  legend.position.inside = c(1, 1),
                  legend.justification = c("right", "top")

            )
        }else{
          psrf_plot <- grid::nullGrob()
        }


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "sigma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "sigma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("sigma"))
  } ## Sigma

  ## Sigma Theta ---------
  if("theta_sd" %in% which.draw){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('theta_sd', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$theta_sd[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('theta_sd', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$theta_sd[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))

      combined_data$var1 <- combined_data$value


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

          median_values <- gelman_plot$shrink[,, "median"]
          upper_values <- gelman_plot$shrink[,, "97.5%"]
          iterations <- gelman_plot$last.iter

          # Create a data frame
          psrf_data <- data.frame(
            Iteration = iterations,
            MedianPSRF = as.vector(median_values),
            UpperPSRF = as.vector(upper_values)
          )

          # Reshape for plotting
          psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
            mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

          psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
            geom_line(aes(linetype = Type, color = Type)) +
            scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
            scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
            labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
            theme_minimal() +
            theme(legend.title = element_blank(),
                  legend.position.inside = c(1, 1),
                  legend.justification = c("right", "top")

            )
        }else{
          psrf_plot <- grid::nullGrob()
        }


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta_sd",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta_sd",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("theta_sd"))
  } ## Sigma theta

  ## z ---------
  if("z" %in% which.draw){

    if(multi_chain){

      draw.item.temp <- apply(draw.item$z, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- nrow(draw.item$z)

      pnames <- c(paste('z [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      nmcmc <- dim(object[[1]]$z)[1]

      for(i in 1:chains){
        z_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$z))
                  for(j in 1:nrow(draw.item$z)){

            iz = draw.item$z[j,1]
            idim = draw.item$z[j,2]

            z_temp[,j] <- object[[i]]$z[,iz,idim]
          }
        chain_list[[i]] <- matrix(z_temp, ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{

      draw.item.temp <- apply(draw.item$z, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- nrow(draw.item$z)

      pnames <- c(paste('z [', draw.item.temp, ']', sep = ''))

      chain_list <- list()
      nmcmc <- dim(object$z)[1]

      z_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$z))

      for(j in 1:nrow(draw.item$z)){

        iz = draw.item$z[j,1]
        idim = draw.item$z[j,2]

        z_temp[,j] <- object$z[,iz,idim]

      }
      chain_list[[1]] <- matrix(z_temp, ncol = length(pnames),
                                dimnames= list(NULL,pnames))


    }

    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:draw.item.num){
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))

      combined_data$var1 <- combined_data$value

      if(multi_chain){

        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position.inside = c(1, 1),
                legend.justification = c("right", "top")
          )

        }else{
          psrf_plot <- grid::nullGrob()
        }

        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))
        if(porder[porder$param == "z",]$idx == max(porder$idx)){
          if(((draw.item.num) > 1)&(i < (draw.item.num))){
            if(interactive()) readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
        
      }else{

        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "z",]$idx == max(porder$idx)){
          if(((draw.item.num) > 1)&(i < (draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }

    which.draw.temp <- setdiff(which.draw.temp, c("z"))
  } ## z

  ## w ---------
  if("w" %in% which.draw){

    if(multi_chain){

      draw.item.temp <- apply(draw.item$w, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- length(draw.item.temp)

      pnames <- c(paste('w [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      nmcmc <- dim(object[[1]]$w)[1]

      for(i in 1:chains){
        w_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$w))
                  for(j in 1:nrow(draw.item$w)){

            iw = draw.item$w[j,1]
            idim = draw.item$w[j,2]

            w_temp[,j] <- object[[i]]$w[,iw,idim]
          }
        chain_list[[i]] <- matrix(w_temp, ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{

      draw.item.temp <- apply(draw.item$w, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- nrow(draw.item$w)

      pnames <- c(paste('w [', draw.item.temp, ']', sep = ''))

      chain_list <- list()
      nmcmc <- dim(object$w)[1]

      w_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$w))

      for(j in 1:nrow(draw.item$w)){

        iw = draw.item$w[j,1]
        idim = draw.item$w[j,2]

        w_temp[,j] <- object$w[,iw,idim]

      }
      chain_list[[1]] <- matrix(w_temp, ncol = length(pnames),
                                dimnames= list(NULL,pnames))


    }

    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:draw.item.num){
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))

      combined_data$var1 <- combined_data$value

      if(multi_chain){

        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)

        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position.inside = c(1, 1),
                legend.justification = c("right", "top")
          )

        }else{
          psrf_plot <- grid::nullGrob()
        }

        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "w",]$idx == max(porder$idx)){
          if(((draw.item.num) > 1)&(i < (draw.item.num))){
            if(interactive()) readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })

        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "w",]$idx == max(porder$idx)){
          if(((draw.item.num) > 1)&(i < (draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }

    which.draw.temp <- setdiff(which.draw.temp, c("w"))
  } ## w
  
  ## zw.dist ---------
  if("zw.dist" %in% which.draw){
    
    if(multi_chain){
      
      
      draw.item.temp <- apply(draw.item$zw.dist, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- length(draw.item.temp)
      
      pnames <- c(paste('zw.dist [', draw.item.temp, ']', sep = ''))
      chains <- object$chains
      
      chain_list <- list()
      nmcmc <- dim(object[[1]]$w)[1]
      
      for(i in 1:chains){
        dist_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$zw.dist))
        for(j in 1:nrow(draw.item$zw.dist)){
          
          iz = draw.item$zw.dist[j,1]
          iw = draw.item$zw.dist[j,2]

          if(dim(object[[i]]$z)[3] == 1){
            dist_temp[,j] <- (object[[i]]$z[,iz,1] - object[[i]]$w[,iw,1])^2
          }else{
            delta <- object[[i]]$z[,iz,] - object[[i]]$w[,iw,]
            dist_temp[,j] <- rowSums(delta^2)
          }
        }
        chain_list[[i]] <- matrix(dist_temp, ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{
      
      draw.item.temp <- apply(draw.item$zw.dist, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- length(draw.item.temp)
      
      pnames <- c(paste('zw.dist [', draw.item.temp, ']', sep = ''))
      
      chain_list <- list()
      nmcmc <- dim(object$w)[1]
      
      dist_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$zw.dist))
      
      for(j in 1:nrow(draw.item$zw.dist)){
        
        iz = draw.item$zw.dist[j,1]
        iw = draw.item$zw.dist[j,2]

        if(dim(object$z)[3] == 1){
          dist_temp[,j] <- (object$z[,iz,1] - object$w[,iw,1])^2
        }else{
          delta <- object$z[,iz,] - object$w[,iw,]
          dist_temp[,j] <- rowSums(delta^2)
        }
        
      }
      chain_list[[1]] <- matrix(dist_temp, ncol = length(pnames),
                                dimnames= list(NULL,pnames))
      
      
    }
    
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]
    
    plots_list <- list()
    
    for(i in 1:draw.item.num){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])
      
      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        y <- as.numeric(chain_data[[j]])
        data.frame(iteration = seq_along(y), value = y, Chain = j)
      }))

      combined_data$var1 <- combined_data$value
      
      
      if(multi_chain){
        
        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")
        
        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")
        
        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })
        
        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        
        if(gelman.diag){
          gelman_plot <- .gelman_plot_safe(mcmclist, param)
        
        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter
        
        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )
        
        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        
        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position.inside = c(1, 1), 
                legend.justification = c("right", "top")  
                
          )
        
        
        }else{
          psrf_plot <- grid::nullGrob()
        }

        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))
        
        if(porder[porder$param == "zw.dist",]$idx == max(porder$idx)){
          if(((draw.item.num) > 1)&(i < (draw.item.num))){
            if(interactive()) readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
        
      }else{
        
        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")
        
        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")
        
        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = as.numeric(acf_data$lag), ACF = as.numeric(acf_data$acf), Chain = j)
        })
        
        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        
        
        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))
        
        
        if(porder[porder$param == "zw.dist",]$idx == max(porder$idx)){
          if(((draw.item.num) > 1)&(i < (draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          if(interactive()) readline(prompt = "Hit <Return> to see next plot")
        }
        
      }
    }
    
    # coda::gelman.diag(mcmclist)
    which.draw.temp <- setdiff(which.draw.temp, c("zw.dist"))
  } ## Dist

  if(gelman.diag == TRUE){
    if(object$chains > 1){
      if(length(chain_list_all) == 0 || any(sapply(chain_list_all, is.null)) || any(sapply(chain_list_all, function(x) length(x) == 0))){
        message("\nWarning: No MCMC samples collected for Gelman-Rubin diagnostic.")
        message("This may occur if no parameters were selected in draw.item.")
      }else{
        cat("\n")
        cat("==================================================\n")
        cat("Gelman and Rubin's convergence diagnostic\n")
        cat("==================================================\n\n")
        mcmc_chains <- lapply(chain_list_all, function(x) coda::mcmc(x, thin = 1))
        mcmclist <- coda::as.mcmc.list(mcmc_chains)
        print(coda::gelman.diag(mcmclist))
        cat("\n")
      }
    }else{
      stop("Gelman and Rubin's convergence diagnostic requires at least two chains for computation.")
    }
  }


}


