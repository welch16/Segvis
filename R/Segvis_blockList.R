#' @title Create a segvis_blockList object
#' @description Constructor for the segvis_blockList class
#' @param ... Object of the segvis_block class
#' @return A segvis_blockList object
#' @export
#' @seealso \code{\link{segvis_block-class}}
Segvis_blockList <- function(...)
{
  aslist = list(...)
  if(length(aslist)==1){
    if(class(aslist[[1]])=="list"){
      segvislist = new("segvis_blockList",aslist[[1]])      
    }else{
      segvislist = new("segvis_blockList",list(aslist[[1]]))
    }    
  }else{
    segvislist = new("segvis_blockList",aslist)
  }
  classes = lapply(segvislist,FUN=class)
  stopifnot(all(classes == "segvis_block"))
  return(segvislist)  
}
