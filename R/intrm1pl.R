#' 1pl LSIRM model with multiplicative effect
#' 
#' @description \link{intrm1pl} integrates functions related to 1pl LSIRM with multiplicative effect. Different missing mechanism can be specified.
#' 
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param missing_data The assumed missing type. One of NA, "mar" and "mcar". Default uses NA. 
#'
#' 
#' @return \code{intrm1pl} returns an object of list. See corresponding function.
#' 
#' @seealso \code{\link{intrm1pl_o}}, \code{\link{intrm1pl_mar}}, \code{\link{intrm1pl_mcar}}
#' @export
intrm1pl = function(data, missing_data = NA, ...) {
  warning("The options \"spikenslab\" and \"fixed_gamma\" are not implemented here because the multiplicative effect does not have distance weight term. It will be ignored if used.")
  if(is.na(missing_data) == TRUE){
    output <- intrm1pl_o(..., as.matrix(data))
    
  }else if(missing_data == 'mar'){
    output <- intrm1pl_mar(..., as.matrix(data))
    
  }else if(missing_data == 'mcar'){
    output <- intrm1pl_mcar(..., as.matrix(data))
  }
  return(output)
}
