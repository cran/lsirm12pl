#' 1pl LSIRM model with normal likelihood
#' 
#' @description \link{lsirm1pl_normal} integrates all functions related to 1pl LSIRM with normal likelihood using multiplicative effect.
#' 
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param spikenslab Whether to use a model selection approach. 
#' @param fixed_gamma Whether fix gamma to 1.
#' @param missing_data The assumed missing type. One of NA, "mar" and "mcar". Default uses NA. 
#'
#' @return \code{lsirm1pl_normal} returns an object of list. See corresponding function.
#' 
#' @note If both \code{spikenslab} and \code{fixed_gamma} are set \code{TRUE}, it returns error because both are related to \code{gamma}.
#' 
#' @seealso \code{\link{lsirm1pl_normal_o}}, \code{\link{lsirm1pl_normal_fixed_gamma}}, \code{\link{lsirm1pl_normal_mar}},
#'
#'  \code{\link{lsirm1pl_normal_mcar}}, \code{\link{lsirm1pl_normal_fixed_gamma_mar}}, \code{\link{lsirm1pl_normal_fixed_gamma_mcar}},
#' 
#'  \code{\link{lsirm1pl_normal_ss}}, \code{\link{lsirm1pl_normal_mar_ss}}, \code{\link{lsirm1pl_normal_mcar_ss}}
#' @export
lsirm1pl_normal = function(data, spikenslab = FALSE, fixed_gamma = FALSE, missing_data = NA, ...) {
  
  if(spikenslab == FALSE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    output <- lsirm1pl_normal_o(..., data)
    
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & is.na(missing_data) == TRUE){
    output <- lsirm1pl_normal_fixed_gamma(..., data)
    
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mar'){
    output <- lsirm1pl_normal_mar(..., data)
    
  }else if(spikenslab == FALSE & fixed_gamma == FALSE & missing_data == 'mcar'){
    output <- lsirm1pl_normal_mcar(..., data)
    
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mar'){
    output <- lsirm1pl_normal_fixed_gamma_mar(..., data)
    
  }else if(spikenslab == FALSE & fixed_gamma == TRUE & missing_data == 'mcar'){
    output <- lsirm1pl_normal_fixed_gamma_mcar(..., data)
    
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & is.na(missing_data) == TRUE){
    output <- lsirm1pl_normal_ss(..., data)
    
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mar'){
    output <- lsirm1pl_normal_mar_ss(..., data)
    
  }else if(spikenslab == TRUE & fixed_gamma == FALSE & missing_data == 'mcar'){
    output <- lsirm1pl_normal_mcar_ss(..., data)
  }else{
    stop('"spikenslab" and "fixed_gamma" cannot both be TRUE.')
  }
  
  return(output)
}
