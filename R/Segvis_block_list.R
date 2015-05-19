#' @title Create a segvis_block_list object
#' @description Constructor for the segvis_block_list class
#' @param ... Object of the segvis_block class
#' @return A segvis_block_list object
#' @export
#' @seealso \code{\link{segvis_block-class}}
Segvis_block_list <- function(...)
{
  aslist <- list(...)
  if(length(aslist) == 1){
    if(class(aslist[[1]]) == "list"){
      segvislist <- new("segvis_block_list",aslist[[1]])      
    }else{
      segvislist <- new("segvis_block_list",list(aslist[[1]]))
    }    
  }else{
    segvislist <- new("segvis_block_list",aslist)
  }
  classes <- lapply(segvislist,FUN=class)
  stopifnot(all(classes == "segvis_block"))
  return(segvislist)  
}
