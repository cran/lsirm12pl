#' 2pl LSIRM model with normal likelihood using multiplicative effect
#' 
#' @description \link{lsirm2pl_normal} integrates all functions related to 2pl LSIRM with normal likelihood.
#' 
#' @param ... Additional arguments for the corresponding function.
#' @param data Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.
#' @param missing_data The assumed missing type. One of NA, "mar", and "mcar". Default uses NA. 
#'
#' @return \code{lsirm2pl_normal} returns an object of list. See corresponding function.
#' 
#' @seealso \code{\link{intrm2pl_normal_o}}, \code{\link{intrm2pl_normal_mar}}, \code{\link{intrm2pl_normal_mcar}}
#' @export
intrm2pl_normal = function(data, missing_data = NA, ...) {
  
  if(is.na(missing_data) == TRUE){
    output <- intrm2pl_normal_o(..., as.matrix(data))
    
  }else if(missing_data == 'mar'){
    output <- intrm2pl_normal_mar(..., as.matrix(data))
    
  }else if(missing_data == 'mcar'){
    output <- intrm2pl_normal_mcar(..., as.matrix(data))
  }
  return(output)
}
