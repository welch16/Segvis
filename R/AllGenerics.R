
# Generic methods for profile class

##  Get methods

#' name
#'
#' @param profile object
#' @return chracter. The name of the profile object
#' @export
#' docType methods
setGeneric("name",
  function(object,...)
  standardGeneric("name")
)           

#' regions
#'
#' @param profile object
#' @return GRangesList. The regions used to calculate the coverage plots, as a GRangesList separated by chromosome
#' @export
#' @docType methods
setGeneric("regions",
  function(object,...)
  standardGeneric("regions")
)

#' bedfiles
#'
#' @param profile object
#' @return chracter. A chracter vector with the names of the bedfiles used to create the readsList object
#' @export
#' docType methods
setGeneric("bedfiles",
  function(object,...)
  standardGeneric("bedfiles")           
)

#' fragLen
#'
#' @param profile object
#' @return numeric. A numeric value representing the width of the extended fragment reads
#' @export
#' docType methods
setGeneric("fragLen",
  function(object,...)
  standardGeneric("fragLen")           
)

#' bandwidth
#'
#' @param profile object
#' @return numeric. A numeric value representing the bandwidth used to smooth the average coverage plot
#' @export
#' docType methods
setGeneric("bandwidth",
  function(object,...)
  standardGeneric("bandwidth")           
)

#' readsList
#'
#' @param profile object
#' @return list. A list made off reads objects. One list for each replicate. Must coincide with size of bedfiles
#' @export
#' docType methods
setGeneric("readsList",
  function(object,...)
  standardGeneric("readsList")           
)

#' matchList
#'
#' @param profile object
#' @return list. A list made off match objects. One list for each replicate. Must coincide with size of bedfiles
#' @export
#' docType methods
#' @rdname profile-methods
setGeneric("matchList",
  function(object,...)
  standardGeneric("matchList")           
)

## Set methods

#' setName
#'
#' @param profile object
#' @param newName character
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setName",
  function(object,newName)
  standardGeneric("setName")           
)

#' setFragLen
#'
#' @param profile object
#' @param newFragLen Numeric value, must be greater of equal to zero
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setFragLen",
  function(object,newFragLen)
  standardGeneric("setFragLen")
)           
                  
#' setBandwidth
#'
#' @param profile object
#' @param newBandwidth Numeric value, must be greater of equal than 1
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setBandwidth",
  function(object,newBandwidth)
  standardGeneric("setBandwidth")
)           
