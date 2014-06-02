#' Sort a GRanges object by the order induced by the start position of the object.
#'
#' @param regions.chr - A GRanges object, for only one chromosome
#' return sorted GRanges object

sort_by_start <- function(gr){

  stopifnot(class(gr) == "GRanges")
  return(gr[sort(start(gr),index.return = TRUE)$ix,])
}
