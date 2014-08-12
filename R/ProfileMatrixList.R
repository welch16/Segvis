#' @title Create a profileMatrixList object
#' @description Constructor for the profileMatrixList class
#' @param ... Objects of the profileMatrix class
#' @return A profileMatrixList object
#' @export
#' @seealso \code{\link{profileMatrixList-class}}
ProfileMatrixList <- function(...){  
  if(class(...) == "list"){
    classes = lapply(...,FUN = class)
    mat = new("profileMatrixList",...)
  }else{
    l = list(...)
    classes = lapply(l,FUN = class)
    mat = new("profileMatrix",l)
  }
  stopifnot(all(classes == "profileMatrix"))
  return(mat)
}
