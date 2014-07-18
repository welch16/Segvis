#' Extend a the reads on a GRanges object to side indicated by it's strand
#'
#' @param gr - A GRanges object
#' @param fragLen - Numeric value of the size to be of the reads on gr
#' @return A Granges object with extended reads
#' @export

extend_to_strand <- function(gr,fragLen)
{
  stopifnot(is(gr,"GRanges"))
  stopifnot(fragLen >=0)
  start(gr) = ifelse(strand(gr) == "-",end(gr)-fragLen +1,start(gr))
  end(gr) = ifelse(strand(gr) == "+",start(gr)+fragLen -1,end(gr))
  return(gr)
}
