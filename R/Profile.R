#' Create a profile object
#'
#' @param regionName - Name of the profile object
#' @param files - Character vector with files that contain the reads
#' @param fileFormat - Character vector with the file format used
#' @param fl - Numeric value of the fragment length used to extend the reads
#' @param maxBw - Numeric value of the maximum possible value to smooth the profile curve
#' @export
#' @return A profile object

Profile <- function(regionName,files,fileF,maxBw,fl)
{
  return(new("profile",name = regionName,files = files,fileFormat = tolower(fileF),
    maxBandwidth = maxBw, fragLen = fl))
}
