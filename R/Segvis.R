#' @title Creates a segvis objec
#'
#' @description Constructor for segvis class
#'
#' @param regionName Name of the segvis object
#'
#' @param file Character string with the file name that contain the reads
#'
#' @param maxBandwidth Numeric value of the maximum possible value to smooth the profile curve
#'
#' @param fragLen Numeric value of the fragment length used to extend the reads
#'
#' @slot isPET logical, Indicates is the reads come from a PET experiments or a SET experiment
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

Segvis <- function(regionName,file,maxBandwidth,chr,fragLen=0,isPET=FALSE)
{
  stopifnot(is.character(regionName))
  stopifnot(is.character(file))
  stopifnot(is.logical(isPET))
  stopifnot(is.character(chr))
  if(isPET){
    warning("The reads are PET, fragLen is set to zero ")
    fragLen =0
  }
  if(chr == "human"){
    chr = paste0("chr",c(1:22,"X","Y"))
  }else if(chr == "mouse"){
    chr = paste0("chr",c(1:19,"X","Y"))
  }else{
    warning("Chromosome",ifelse(length(chr)==1,"",
                                "s")," supplied by user")
  }  
  return(new("segvis",name = regionName,file=file,
    maxBandwidth = maxBandwidth, fragLen = fragLen,isPET=isPET,chr =chr))
}
