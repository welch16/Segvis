#' @title Create a profileMatrixList object
#' @description Constructor for the profileMatrixList class
#' @param ... Objects of the profileMatrix class
#' @return A profileMatrixList object
#' @export
#' @seealso \code{\link{profileMatrixList-class}}
ProfileMatrixList <- function(...){
  aslist = list(...)
  if(length(aslist)==1){
    if(class(aslist[[1]])=="list"){
      mat = new("profileMatrixList",aslist[[1]])      
    }else{
      mat = new("profileMatrixList",list(aslist[[1]]))
    }    
  }else{
    mat = new("profileMatrixList",aslist)
  }
  classes = lapply(mat,FUN=class)
  stopifnot(all(classes == "profileMatrix"))
  return(mat)
}
