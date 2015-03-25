#' @title Create a segvisList object
#' @description Constructor for the segvisList class
#' @param ... Object of the segvis class
#' @return A segvisList object
#' @export
#' @seealso \code{\link{segvis-class}}
SegvisList <- function(...)
{
  aslist = list(...)
  if(length(aslist)==1){
    if(class(aslist[[1]])=="list"){
      segvislist = new("segvisList",aslist[[1]])      
    }else{
      segvislist = new("segvisList",list(aslist[[1]]))
    }    
  }else{
    segvislist = new("segvisList",aslist)
  }
  classes = lapply(segvislist,FUN=class)
  stopifnot(all(classes == "segvis"))
  return(segvislist)  
}
