#' Print the summary the result of LSIRM
#'
#' @description \link{print.summary.lsirm} is used to print summary the result of LSIRM.
#'
#' @param x List; summary of LSIRM with \code{summary.lsirm}.
#' @param items Numeric or Character vector; the items to display in the parameter summaries. Either a vector of item indices (e.g. \code{1:5}) or item names (e.g. \code{c("item 1", "item 3")}). Default is \code{NULL} (displays all items).
#' @param respondents Numeric or Character vector; the respondents to display in the respondent summaries (theta). Either a vector of respondent indices (e.g. \code{1:10}) or respondent names (e.g. \code{c("respondent 1", "respondent 2")}). Default is \code{NULL} (displays all respondents if theta is requested).
#' @param which Character vector; specifies which parameters to print in the console. Options are \code{"beta"} (item difficulty/threshold), \code{"alpha"} (item discrimination), and \code{"theta"} (respondent ability). Or set to \code{"all"} to print all. Default is \code{NULL} (dynamically prints \code{"beta"} for 1PL models, and both \code{"beta"} and \code{"alpha"} for 2PL models).
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
print.summary.lsirm <- function(x, items = NULL, respondents = NULL, which = NULL, ...){
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
        v_num <- as.numeric(v)
        formatted <- format(round(v_num, digits = digits), nsmall = digits, trim = TRUE)
        formatted[is.na(v_num)] <- "NA"
        return(formatted)
      }

      cat(sprintf("Adapt interval: %s\tAdapt rate: %s\n",
                  x$tuning$adapt_interval, fmt(x$tuning$adapt_rate, 3)))

      params <- c("beta", "theta", "z", "w", "gamma", "alpha")
      
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

  if (is.null(items)) {
    items <- x$items
  }
  if (is.null(respondents)) {
    respondents <- x$respondents
  }
  if (is.null(which)) {
    which <- x$which
  }

  is_ordinal <- grepl("ordinal", x$method, ignore.case = TRUE)
  beta_label <- if (is_ordinal) "beta (item threshold)" else "beta (item difficulty)"

  if (is.null(which)) {
    if (!is.null(x$alpha.coef)) {
      which <- c("beta", "alpha")
    } else {
      which <- c("beta")
    }
  } else if (length(which) == 1 && which == "all") {
    which <- c("beta", "alpha", "theta")
  }

  # 1. Print beta (item difficulty/threshold) if requested
  if ("beta" %in% which) {
    coef_to_print <- x$coef
    subsetted_items <- FALSE
    total_items <- 0

    if (is_ordinal) {
      row_items <- sapply(strsplit(rownames(x$coef), ":"), `[`, 1)
      unique_items <- unique(row_items)
      total_items <- length(unique_items)
      
      if (!is.null(items)) {
        if (is.numeric(items)) {
          valid_indices <- items[items >= 1 & items <= total_items]
          selected_items <- unique_items[valid_indices]
        } else {
          selected_items <- items[items %in% unique_items]
        }
        selected_rows <- row_items %in% selected_items
        coef_to_print <- x$coef[selected_rows, , drop = FALSE]
        subsetted_items <- TRUE
      }
    } else {
      total_items <- nrow(x$coef)
      if (!is.null(items)) {
        if (is.numeric(items)) {
          valid_indices <- items[items >= 1 & items <= total_items]
          coef_to_print <- x$coef[valid_indices, , drop = FALSE]
        } else {
          coef_to_print <- x$coef[rownames(x$coef) %in% items, , drop = FALSE]
        }
        subsetted_items <- TRUE
      }
    }

    if(x$n.chains == 1){
      cat(sprintf("Posterior estimates of %s parameters: \n\n", beta_label))
    }else{
      cat(sprintf("Posterior estimates of %s parameters of chain %d: \n\n", beta_label, x$chain))
    }
    printCoefmat(coef_to_print)

    if (subsetted_items) {
      cat(sprintf("\n[Showing selected items. Showing %d items out of total %d items.]\n", 
                  if (is_ordinal) length(unique(sapply(strsplit(rownames(coef_to_print), ":"), `[`, 1))) else nrow(coef_to_print), 
                  total_items))
    }
    cat("\n")
  }

  # 2. Print alpha (item discrimination) if requested and exists
  if ("alpha" %in% which && is.null(x$alpha.coef)) {
    warning("The model does not contain 'alpha' (item discrimination) parameters. Skipped printing 'alpha'.", call. = FALSE)
  }

  if ("alpha" %in% which && !is.null(x$alpha.coef)) {
    alpha_to_print <- x$alpha.coef
    subsetted_alpha <- FALSE
    total_items <- nrow(x$alpha.coef)

    if (!is.null(items)) {
      if (is.numeric(items)) {
        valid_indices <- items[items >= 1 & items <= total_items]
        alpha_to_print <- x$alpha.coef[valid_indices, , drop = FALSE]
      } else {
        alpha_to_print <- x$alpha.coef[rownames(x$alpha.coef) %in% items, , drop = FALSE]
      }
      subsetted_alpha <- TRUE
    }

    if(x$n.chains == 1){
      cat("Posterior estimates of alpha (item discrimination) parameters: ","\n\n")
    }else{
      cat("Posterior estimates of alpha (item discrimination) parameters of chain", x$chain,": ","\n\n")
    }
    printCoefmat(alpha_to_print)

    if (subsetted_alpha) {
      cat(sprintf("\n[Showing selected items. Showing %d items out of total %d items.]\n", 
                  nrow(alpha_to_print), total_items))
    }
    cat("\n")
  }

  # 3. Print theta (respondent ability) if requested
  if (("theta" %in% which || "respondent" %in% which) && !is.null(x$theta.coef)) {
    theta_to_print <- x$theta.coef
    subsetted_resp <- FALSE
    total_respondents <- nrow(x$theta.coef)

    if (!is.null(respondents)) {
      if (is.numeric(respondents)) {
        valid_resp_indices <- respondents[respondents >= 1 & respondents <= total_respondents]
        theta_to_print <- x$theta.coef[valid_resp_indices, , drop = FALSE]
      } else {
        theta_to_print <- x$theta.coef[rownames(x$theta.coef) %in% respondents, , drop = FALSE]
      }
      subsetted_resp <- TRUE
    }

    if(x$n.chains == 1){
      cat("Posterior estimates of theta (respondent ability) parameters: ","\n\n")
    }else{
      cat("Posterior estimates of theta (respondent ability) parameters of chain", x$chain,": ","\n\n")
    }
    printCoefmat(theta_to_print)

    if (subsetted_resp) {
      cat(sprintf("\n[Showing selected respondents. Showing %d respondents out of total %d respondents.]\n", 
                  nrow(theta_to_print), total_respondents))
    }
    cat("\n")
  }

  cat("\n---------------------------","\n\n")
  cat("Overall BIC (Smaller is better) :",x$BIC,"\n")
  cat("\nMaximum Log-posterior Iteration: ", "\n")
  printCoefmat(x$map.inf)
}
