#' Print the summary the result of LSIRM
#'
#' @description \link{print.summary.lsirm} is used to print summary the result of LSIRM.
#'
#' @param x List; summary of LSIRM with \code{summary.lsirm}.
#' @param ... Additional arguments.
#'
#' @return \code{print.summary.lsirm} return a summary of LSIRM.
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl())
#' summary(lsirm_result)
#' }
#' @rdname print.summary.lsirm
#' @export
print.summary.lsirm <- function(x, ...){
  cat("==========================","\n")
  cat("Summary of model","\n")
  cat("==========================","\n\n")
  cat("Call:\t", deparse(x$call), "\n")
  cat("Model:\t", x$method, "\n")
  cat("Data type:\t", x$dtype, "\n")
  cat("Variable Selection:\t", x$ss, "\n")
  cat("Missing:\t", x$missing, "\n")
  cat(sprintf("MCMC sample of size %i, after burnin of %i iteration",
              x$mcmc.opt$niter,x$mcmc.opt$nburn),"\n\n")

  if(!is.null(x$tuning)){
    use_adapt <- isTRUE(x$tuning$use_adapt)
    cat("Adaptive MCMC (burn-in only):\t", use_adapt, "\n")
    if(use_adapt){
      fmt <- function(v, digits = 3){
        if(length(v) == 0 || is.null(v)) return("NA")
        # Ensure input is numeric; NAs will propagate
        v_num <- as.numeric(v)
        # round and format; handling vector input
        formatted <- format(round(v_num, digits = digits), nsmall = digits, trim = TRUE)
        formatted[is.na(v_num)] <- "NA"
        return(formatted)
      }

      cat(sprintf("Adapt interval: %s\tAdapt rate: %s\n",
                  x$tuning$adapt_interval, fmt(x$tuning$adapt_rate, 3)))

      # Create a summary table for adaptive MCMC
      params <- c("beta", "theta", "z", "w", "gamma", "alpha")
      
      # Extract values with helper to handle various NA/NULL types
      get_val <- function(key, default=NA) {
        if(is.null(x$tuning[[key]])) return(default)
        if(length(x$tuning[[key]]) == 0) return(default)
        val <- x$tuning[[key]]
        if(is.na(val)) return(default)
        return(as.numeric(val))
      }

      targets <- c(get_val("target_accept_beta"), get_val("target_accept_theta"), 
                   get_val("target_accept_zw"), get_val("target_accept_zw"), 
                   get_val("target_accept_gamma"), get_val("target_accept_alpha"))
      
      jump_init <- c(get_val("jump_beta_init"), get_val("jump_theta_init"),
                     get_val("jump_z_init"), get_val("jump_w_init"),
                     get_val("jump_gamma_init"), get_val("jump_alpha_init"))

      jump_final <- c(get_val("jump_beta_final"), get_val("jump_theta_final"),
                      get_val("jump_z_final"), get_val("jump_w_final"),
                      get_val("jump_gamma_final"), get_val("jump_alpha_final"))
      
      acc_burn <- c(get_val("accept_beta_burn"), get_val("accept_theta_burn"),
                    get_val("accept_z_burn"), get_val("accept_w_burn"),
                    get_val("accept_gamma_burn"), get_val("accept_alpha_burn"))

      acc_last <- c(get_val("accept_beta_lastwin"), get_val("accept_theta_lastwin"),
                    get_val("accept_z_lastwin"), get_val("accept_w_lastwin"),
                    get_val("accept_gamma_lastwin"), get_val("accept_alpha_lastwin"))

      # Filter out rows where all values are NA (e.g. alpha in 1PL)
      # Check if a parameter is relevant: at least one of target, jump, or accept is not NA
      is_relevant <- !is.na(targets) | !is.na(jump_init) | !is.na(jump_final) | !is.na(acc_burn) | !is.na(acc_last)
      
      if(any(is_relevant)) {
        adapt_df <- data.frame(
          Parameter = params[is_relevant],
          Target = fmt(targets[is_relevant]),
          Jump_Init = fmt(jump_init[is_relevant]),
          Jump_Final = fmt(jump_final[is_relevant]),
          Acc_BurnIn = fmt(acc_burn[is_relevant]),
          Acc_LastWin = fmt(acc_last[is_relevant])
        )
        print(adapt_df, row.names = FALSE)
      }

    }
    cat("\n")
  }
  if(x$n.chains == 1){
    cat("Posterior estimates of item parameters: ","\n\n")
  }else{
    cat("Posterior estimates of item parameters of chain", x$chain,": ","\n\n")
  }
  printCoefmat(x$coef)
  cat("\n---------------------------","\n\n")
  cat("Overall BIC (Smaller is better) :",x$BIC,"\n")
  cat("\nMaximum Log-posterior Iteration: ", "\n")
  printCoefmat(x$map.inf)
}
