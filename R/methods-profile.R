
# Methods for profile class

## Get methods

#' name
#'
#' @param profile object
#' @return character. The name of the object
#' @docType methods
#' @rdname profile-methods
setMethods("name",
  signature = signature(object = "profile"),
  definition = function(object)object@name
)           

#' regions
#'
#' @param profile object
#' @return GRangesList. The regions used to calculate the coverage plots, as a GRangesList separated by chromosome
#' @docType methods
#' @rdname profile-methods
setMethods("regions",
  signature = signature(object = "profile"),
  definition = function(object)object@regions
)           

#' bedfiles
#'
#' @param profile object
#' @return chracter. A chracter vector with the names of the bedfiles used to create the readsList object
#' @docType methods
#' @rdname profile-methods
setMethods("bedfiles",
  signature = signature(object ="profile"),
  definition = function(object)object@bedfiles
)           

#' fragLen
#'
#' @param profile object
#' @return numeric. A numeric value representing the width of the extended fragment reads
#' @docType methods
#' @rdname profile-methods
setMethods("fragLen",
  signature = signature(object = "profile"),
  definition = function(object)object@fragLen
)           


#' bandwidth
#'
#' @param profile object
#' @return numeric. A numeric value representing the bandwidth used to smooth the average coverage plot
#' @docType methods
#' @rdname profile-methods
setMethods("bandwidth",
  signature = signature(object = "profile"),
  definition = function(object)object@bandwidth
)           

#' readsList
#'
#' @param profile object
#' @return list. A list made off reads objects. One list for each replicate. Must coincide with size of bedfiles
#' docType methods
#' @rdname profile-methods
setMethods("readsList",
  signature = signature(object = "profile"),
  definition = function(object)object@readsList           
)
           
#' matchList
#'
#' @param profile object
#' @return list. A list made off match objects. One list for each replicate. Must coincide with size of bedfiles
#' @docType methods
#' @rdname profile-methods
setMethods("matchList",
  signature = signature(object = "profile"),
  definition = function(object)object@matchList
)

## Set methods

#' setName
#'
#' @param profile object
#' @param newName character
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setName",
  signature = signature(object = "profile",newName = "character"),
  definition = function(object,newName){
    object@name = newName
    return(object)
})

#' setRegions
#'
#' @param profile object
#' @param newRegions GRangesList
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setRegions",
  signature = signature(object = "profile",newRegions = "GRangesList"),
  definition = function(object,newRegions){
    stopifnot(class(newRegions) == "GRangesList")
    object@regions = newRegions
    return(object)
})

#' setFragLen
#'
#' @param profile object
#' @param newFragLen Numeric value, must be greater of equal to zero
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setFragLen",
  signature = signature(object = "profile",newFragLen = "numeric"),
  definition = function(object,newFragLen){
    stopifnot(newFragLen >= 0)
    stopifnot(newFragLen == floor(newFragLen))
    object@fragLen = newFraLen
    return(object)
})    

#' setBandwidth
#'
#' @param profile object
#' @param newBandwidth Numeric value, must be greater of equal to one
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setBandwidth",
  signature = signature(object = "profile",newBandwidth = "numeric"),
  definition = function(object,newBandwidth){
    stopifnot(newBandwidth >= 1)
    stopifnot(newBandwidth == floor(newBandwidth))
    object@bandwidth = newBandwidth
    return(object)
})    

#' show
#'
#' @param profile object
#' @docType methods
#' @rdname profile-methods
setMethods("show",
  signature = signature(object = "profile"),
  definition = function(object){
    cat("---------------------------\n")
    cat("Profile for",name(object),"peaks\n")
    cat("Fragment length:",fragLen(object),"\n")
    cat("Bandwidth:",bandwidth(object),"\n")
    if(q <- length(regions(object)) > 0){
      cat("Using regions for",length(regions(object)),"chromosomes\n")
    }else{
      cat("**Not regions loaded**\n")
    }     
    cat("Using reads files:\n")
    cat(bedfiles(object),sep = "\n")
    cat("---------------------------\n")
})
          
