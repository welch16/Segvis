#' @title Create a profileMatrixList object
#' @description Constructor for the profileMatrixList class
#' @param ... Objects of the profileMatrix class
#' @return A profileMatrixList object
#' @export
#' @seealso \code{\link{profileMatrixList-class}}
ProfileMatrixList <- function(...){
  l = list(...)
  classes = lapply(l,FUN = class)
  stopifnot(all(classes == "profileMatrix"))
  return(new("profileMatrixList",l))
}
