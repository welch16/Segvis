#' @title Create a profileMatrix object
#' @description Constructor for profileMatrix class
#' @param name Character with the name of the profileMatrix object
#' @param regions GRanges object with the regions considered for the profile, the width of all regions must be the same
#' @export
#' @return a profileMatrix object

ProfileMatrix <- function(name,regions)
{
  stopifnot(class(name) == "character")
  stopifnot(class(regions) == "GRanges")
  return(new("profileMatrix",name = name,regions = regions))
}
