#' @title Creates a profile objec
#' @description Constructor for profile class
#' @param regionName Name of the profile object
#' @param file Character string with the file name that contain the reads
#' @param fileFormat Character vector with the file format used
#' @param fragLen Numeric value of the fragment length used to extend the reads
#' @param maxBandwidth Numeric value of the maximum possible value to smooth the profile curve
#' @param remChr Character Vector with the chromosomes to be ignored
#' @export
#' @return A profile object
#' @examples
#' rn = "sites"
#' f = "somefile.bam";fileF = "bam"
#' maxBw = 1;fl = 200
#' Profile(regionName = rn,file =f,fileFormat=fileF,maxBandwidth = maxBw, fragLen = fl)

Profile <- function(regionName,file,fileFormat,maxBandwidth,fragLen,remChr=NULL)
{
  stopifnot(is.character(regionName))
  stopifnot(is.character(file))
  stopifnot(is.character(fileFormat))
  if(is.null(remChr))remChr=""
  stopifnot(is.character(remChr))
  return(new("profile",name = regionName,file=file,fileFormat = tolower(fileFormat),
    maxBandwidth = maxBandwidth, fragLen = fragLen,remChr = remChr))
}
