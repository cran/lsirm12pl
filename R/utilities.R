
check.datatype <- function(data, missing.val = 99, ...){
  if(exists("missing.val")){ # If missing.var is imputed
    missing.temp = missing.val
  }else{missing.temp = 99}

  check.v <- c(0, 1, missing.temp)
  check.r <- length(setdiff(unique(c(data.matrix(data))), check.v)) == 0 # If the unique values are exist except check.v
  return(check.r) # TRUE: Binary, FALSE: Continuous

  # if(!exists('missing_data')){
  #   missing.temp = c()
  # }else{
  #   missing.temp = missing.var
  # }
  # check.v <- c(0,1,missing.temp)
  # check.r <- sum(check.v %in% unique(c(data.matrix(data)))) == length(check.v)  # If data have (0,1,90), cannot catch continuous
  # return(check.r)
}

