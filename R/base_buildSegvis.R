#' Builds a segvis object
#'
#' Constructor for segvis class
#'
#' @param name Name of the segvis object
#'
#' @param file Character string with the file name that contain the reads
#'
#' @param maxBandwidth Numeric value of the maximum possible value to smooth the profile curve
#'
#' @param fragLen Numeric value of the fragment length used to extend the reads
#'
#' @param isPET logical, Indicates is the reads come from a PET experiments or a SET experiment
#'
#' @export
#'
#' @return A segvis object
#'
#' @rdname buildSegvis
#' @name buildSegvis
#' @examples
#' \dontrun{
#' rn = "sites"
#' f = "somefile.bam"
#' maxBw = 1;fl = 200
#' buildSegvis(regionName = rn,file =f,maxBandwidth = maxBw, fragLen = fl)
#' }

buildSegvis <- function(name,file,maxBandwidth,chr,fragLen=0,isPET=FALSE)
{
  stopifnot(is.character(name))
  stopifnot(is.character(file))
  stopifnot(is.logical(isPET))
  stopifnot(is.character(chr))
  if(isPET){
    warning("The reads are PET, fragLen is set to zero ")
    fragLen <- 0
  }
  if(length(chr) == 1){
    if(chr == "human"){
      chr <- paste0("chr",c(1:22,"X","Y"))
    }else if(chr == "mouse"){
      chr <- paste0("chr",c(1:19,"X","Y"))
    }else{
      warning("Chromosome supplied by user")
    }    
  }else{
    warning("Chromosomes supplied by user")
  }
  return(new("segvis",name = name,file=file,
    maxBandwidth = maxBandwidth, fragLen = fragLen,isPET=isPET,chr =chr))
}
