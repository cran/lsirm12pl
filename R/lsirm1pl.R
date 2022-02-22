#' 1pl LSIRM model
#' 
#' @description \link{lsirm1pl} integrates all functions related to 1pl LSIRM
#' 
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param spikenslab Whether to use a model selection approach. 
#' @param fixed_gamma Whether fix gamma to 1.
#' @param missing_data The assumed missing type. One of NA, "mar" and "mcar". Default uses NA. 
#'
#' @return \code{lsirm1pl} returns an object of list. See corresponding function.
#' 
#' @note If both \code{spikenslab} and \code{fixed_gamma} are set \code{TRUE}, it returns error because both are related to \code{gamma}.
#' 
#' @seealso \code{\link{lsirm1pl_o}}, \code{\link{lsirm1pl_fixed_gamma}}, \code{\link{lsirm1pl_mar}},  \code{\link{lsirm1pl_mcar}},\code{\link{lsirm1pl_fixed_gamma_mar}}, \code{\link{lsirm1pl_fixed_gamma_mcar}}, \code{\link{lsirm1pl_ss}}, \code{\link{lsirm1pl_mar_ss}}, \code{\link{lsirm1pl_mcar_ss}}
#' 
#' @export
lsirm1pl = function(data, spikenslab = FALSE, fixed_gamma = FALSE, missing_data = NA, ...) {
  
  if(spikenslab == FALSE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    output <- lsirm1pl_o(..., as.matrix(data))
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mar'){
    output <- lsirm1pl_mar(..., as.matrix(data))
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mcar'){
    output <- lsirm1pl_mcar(..., as.matrix(data))
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & is.na(missing_data) == TRUE){
    output <- lsirm1pl_fixed_gamma(..., as.matrix(data) )
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mar'){
    output <- lsirm1pl_fixed_gamma_mar(..., as.matrix(data))
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mcar'){
    output <- lsirm1pl_fixed_gamma_mcar(..., as.matrix(data))
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    output <- lsirm1pl_ss(..., as.matrix(data))
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mar'){
    output <- lsirm1pl_mar_ss(..., as.matrix(data) )
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mcar'){
    output <- lsirm1pl_mcar_ss(..., as.matrix(data))
  }else{
    stop('The options "spikenslab" and "fixed_gamma" cannot be set TRUE at the same time.')
  }
  
  
  return(output)
}
