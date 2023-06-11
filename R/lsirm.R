#' Fit a LSIRM model
#'
#' @description \link{lsirm} is used to fit 1PL LSIRM and 2PL LSIRM using Bayesian method as described in Jeon et al. (2021).
#'
#' @param data Matrix; binary or continuous item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param ... Additional arguments for the corresponding function.
#'
#' @return \code{lsirm} returns an object of list. See corresponding function.
#' @seealso \link{lsirm1pl} and \link{lsirm2pl}
#' @export
lsirm = function(data, ...) UseMethod("lsirm")

#' Formula function for LSIRM model
#' @description \link{lsirm.formula} is formula object.
#'
#' @param formula The form of formula is \code{lsirm(A ~ <term 1>(<term 2> + <term 3> ...))}, where \code{A} is binary or continuous item response matrix to be analyzed, \code{<term1>} is the model you want to fit and has one of the following values: "lsirm1pl" and "lsirm2pl"., and \code{<term 2>}, \code{<term 3>}, etc., are each option for the model.
#' @param ... Additional arguments for the corresponding function.
#'
#' @export
lsirm.formula = function(formula, ...){
  data <- get(eval(formula)[[2]])
  argument <- rlang::call_args(eval(formula)[[3]])
  argument$data <- data
  func <- as.character(eval(formula)[[3]][[1]])
  output <- do.call(func,argument)
  output$call <- match.call()
  class(output) <- "lsirm"
  return(output)
}

#' Fit a 1pl LSIRM model
#'
#' @description \link{lsirm1pl} integrates all functions related to 1pl LSIRM
#'
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary or continuous item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param spikenslab Whether to use a model selection approach.
#' @param fixed_gamma Whether fix gamma to 1.
#' @param missing_data The assumed missing type. One of NA, "mar" and "mcar". Default uses NA.
#'
#' @return \code{lsirm1pl} returns an object of list. See corresponding function.
#'
#' @note If both \code{spikenslab} and \code{fixed_gamma} are set \code{TRUE}, it returns error because both are related to \code{gamma}.
#'
#' @seealso \code{\link{lsirm1pl_o}}, \code{\link{lsirm1pl_fixed_gamma}}, \code{\link{lsirm1pl_mar}},  \code{\link{lsirm1pl_mcar}},\code{\link{lsirm1pl_fixed_gamma_mar}}, \code{\link{lsirm1pl_fixed_gamma_mcar}}, \code{\link{lsirm1pl_ss}}, \code{\link{lsirm1pl_mar_ss}}, \code{\link{lsirm1pl_mcar_ss}}, \code{\link{lsirm1pl_normal_o}}, \code{\link{lsirm1pl_normal_fixed_gamma}}, \code{\link{lsirm1pl_normal_mar}},  \code{\link{lsirm1pl_normal_mcar}},\code{\link{lsirm1pl_normal_fixed_gamma_mar}}, \code{\link{lsirm1pl_normal_fixed_gamma_mcar}}, \code{\link{lsirm1pl_normal_ss}}, \code{\link{lsirm1pl_normal_mar_ss}}, \code{\link{lsirm1pl_normal_mcar_ss}}
#' @export
lsirm1pl = function(data, spikenslab = FALSE, fixed_gamma = FALSE, missing_data = NA, ...) {

  if(spikenslab == FALSE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){ # if missing.val imputed
      output <- lsirm1pl_o(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_o(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_mar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_mar(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_mcar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_mcar(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_fixed_gamma(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_fixed_gamma(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_fixed_gamma_mar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_fixed_gamma_mar(..., data)
      output$dtype <- "continuous"
    }

  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_fixed_gamma_mcar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_fixed_gamma_mcar(..., data)
      output$dtype <- "continuous"
    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_ss(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_ss(..., data)
      output$dtype <- "continuous"
    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_mar_ss(..., data )
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_mar_ss(..., data )
      output$dtype <- "continuous"
    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      output <- lsirm1pl_mcar_ss(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm1pl_normal_mcar_ss(..., data)
      output$dtype <- "continuous"
    }

  }else{
    stop('The options "spikenslab" and "fixed_gamma" cannot be set TRUE at the same time.')
  }
  output$call <- match.call()
  output$method <- "lsirm1pl"
  output$missing <- missing_data
  output$varselect <- spikenslab
  class(output) <- "lsirm"
  return(output)
}



#' Fit a 2pl LSIRM model
#'
#' @description \link{lsirm2pl} integrates all functions related to 2pl LSIRM
#'
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary or continuous item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param spikenslab Whether to use a model selection approach.
#' @param fixed_gamma Whether fix gamma to 1.
#' @param missing_data The assumed missing type. One of NA, "mar" and "mcar". Default uses NA.
#'
#' @return \code{lsirm2pl} returns an object of list. See corresponding function.
#'
#' @note If both \code{spikenslab} and \code{fixed_gamma} are set \code{TRUE}, it returns error because both are related to \code{gamma}.
#'
#' @seealso \code{\link{lsirm2pl_o}}, \code{\link{lsirm2pl_fixed_gamma}}, \code{\link{lsirm2pl_mar}},  \code{\link{lsirm2pl_mcar}},\code{\link{lsirm2pl_fixed_gamma_mar}}, \code{\link{lsirm2pl_fixed_gamma_mcar}}, \code{\link{lsirm2pl_ss}}, \code{\link{lsirm2pl_mar_ss}}, \code{\link{lsirm2pl_mcar_ss}}, \code{\link{lsirm2pl_normal_o}}, \code{\link{lsirm2pl_normal_fixed_gamma}}, \code{\link{lsirm2pl_normal_mar}},  \code{\link{lsirm2pl_mcar}},\code{\link{lsirm2pl_normal_fixed_gamma_mar}}, \code{\link{lsirm2pl_normal_fixed_gamma_mcar}}, \code{\link{lsirm2pl_normal_ss}}, \code{\link{lsirm2pl_normal_mar_ss}}, \code{\link{lsirm2pl_normal_mcar_ss}}
#'
#' @export
lsirm2pl = function(data, spikenslab = FALSE, fixed_gamma = FALSE, missing_data = NA, ...) {

  if(spikenslab == FALSE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_o(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_o(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_mar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_mar(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_mcar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_mcar(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_fixed_gamma(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_fixed_gamma(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_fixed_gamma_mar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_fixed_gamma_mar(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_fixed_gamma_mcar(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_fixed_gamma_mcar(..., data)
      output$dtype <- "continuous"
    }
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_ss(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_ss(..., data)
      output$dtype <- "continuous"
    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mar'){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_mar_ss(..., data )
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_mar_ss(..., data )
      output$dtype <- "continuous"
    }

  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mcar'){
    if(check.datatype(data, ...)){
      output <- lsirm2pl_mcar_ss(..., data)
      output$dtype <- "binary"
    } else{
      output <- lsirm2pl_normal_mcar_ss(..., data)
      output$dtype <- "continuous"
    }

  }else{
    stop('The options "spikenslab" and "fixed_gamma" cannot be set TRUE at the same time.')
  }
  output$call <- match.call()
  output$method <- "lsirm2pl"
  output$missing <- missing_data
  output$varselect <- spikenslab
  class(output) <- "lsirm"
  return(output)
}



