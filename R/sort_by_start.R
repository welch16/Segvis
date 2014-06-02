#' Sort a GRanges object by the order induced by the start position of the object.
#'
#' @param regions.chr - A GRanges object, for only one chromosome
#' return sorted GRanges object

sort_by_start <- function(regions.chr){

  stopifnot(class(regions.chr) == "GRanges")
  return(regions.chr[sort(start(regions.chr),index.return = TRUE)$ix,])
}
