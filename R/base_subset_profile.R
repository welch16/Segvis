#' Programming function to obtain a logical vector of the peaks in regions(object) such that satisfy condition_all
#' @param object profileMatrix object
#' @param condition_call name object with the filtering condition
#' @return A logical vector with the same length as the number of peaks
#' @rdname filter_profileMat
.subset_profileMat_logical<- function(object,condition_call){  
  cond = as.logical(eval(condition_call,as(regions(object),"data.frame"),parent.frame()))
  return(cond)
}

#' Programming function to filter a profileMatrix object
#' @param cond Logical Vector with the peaks that satisfy certain condition
#' @return A profileMat object with the peaks that such that cond is TRUE
#' @rdname filter_profileMat
.filter_profileMat <- function(object,cond){
  mat = profileMat(object)[cond,]
  if(sum(cond)==1){
    mat = t(as.matrix(mat))
  }
  out = new("profileMatrix",name = name(object),regions = regions(object)[cond],
      bandwidth = bandwidth(object),normConst = normConst(object),profileMat = mat,
               .isScaled = object@.isScaled)
  return(out)    
}
