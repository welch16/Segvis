#' Create a profile object
#'
#' @param regionName - Name of the profile object
#' @param files - Character vector with files that contain the reads
#' @param fl - Numeric value of the fragment length used to extend the reads
#' @param bw - Numeric value of the bandwidth used to smooth the coverage
#' @export
#' @return A profile object

create_profile <- function(regionName,files,fl,bw)
{
  return(new("profile",name = regionName,bedfiles = files,
    fragLen = fl, bandwidth = bw))             
}
