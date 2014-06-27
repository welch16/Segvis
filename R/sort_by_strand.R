#' Sort the reads on a GRanges object depending on it's strand
#'
#' @param gr - A GRanges object
#' @param stra - Unique strand of the object
#' @return A GRanges object
#' @export

sort_by_strand <- function(gr,stra)
{
  if(stra == "-"){
    gr = gr[sort(end(gr),index.return = TRUE)$ix]
  }else{
    gr = gr[sort(start(gr),index.return = TRUE)$ix]
  }
  return(gr)
}
