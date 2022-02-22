#' 2pl LSIRM model using multiplicative effect
#' 
#' @description \link{intrm2pl} integrates all functions related to 2pl LSIRM using multiplicative effect.
#' 
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param missing_data The assumed missing type. One of NA, "mar" and "mcar". Default uses NA. 
#'
#' @return \code{intrm2pl} returns an object of list. See corresponding function.
#' 
#' @seealso \code{\link{intrm2pl_o}}, \code{\link{intrm2pl_mar}}, \code{\link{intrm2pl_mcar}}
#' @export
intrm2pl = function(data, missing_data = NA, ...){
  warning("The options \"spikenslab\" and \"fixed_gamma\" are not implemented here because the multiplicative effect does not have distance weight term. It will be ignored if used")
  if(is.na(missing_data) == TRUE){
    output <- intrm2pl_o(..., as.matrix(data))
    
  }else if(missing_data == 'mar'){
    output <- intrm2pl_mar(..., as.matrix(data))
    
  }else if(missing_data == 'mcar'){
    output <- intrm2pl_mcar(..., as.matrix(data))
  }
  return(output)
}
