#' @title Creates a segvis objec
#'
#' @description Constructor for segvis class
#'
#' @param regionName Name of the segvis object
#'
#' @param file Character string with the file name that contain the reads
#'
#' @param fragLen Numeric value of the fragment length used to extend the reads
#'
#' @param maxBandwidth Numeric value of the maximum possible value to smooth the profile curve
#'
#' @export
#'
#' @return A segvis object
#'
#' @examples
#' rn = "sites"
#' f = "somefile.bam"
#' maxBw = 1;fl = 200
#' Segvis(regionName = rn,file =f,maxBandwidth = maxBw, fragLen = fl)

Segvis <- function(regionName,file,fileFormat,maxBandwidth,fragLen)
{
  stopifnot(is.character(regionName))
  stopifnot(is.character(file))
  return(new("segvis",name = regionName,file=file,
    maxBandwidth = maxBandwidth, fragLen = fragLen))
}
