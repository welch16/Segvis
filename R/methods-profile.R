

# Generic methods for profile class

#' regions
#'
#' @param profile object
#' @return GRangesList. The regions used to calculate the coverage plots, as a GRangesList separated by chromosome
#' @export
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
#' @export
#' docType methods
#' @rdname profile-methods
setMethods("bedfiles",
  signature = signature(object ="profile"),
  definition = function(object)object@bedfiles
)           

#' fragLen
#'
#' @param profile object
#' @return numeric. A numeric value representing the width of the extended fragment reads
#' @export
#' docType methods
#' @rdname profile-methods
setMethods("fragLen",
  signature = signature(object = "profile"),
  definition = function(object)object@fragLen
)           


#' bandwidth
#'
#' @param profile object
#' @return numeric. A numeric value representing the bandwidth used to smooth the average coverage plot
#' @export
#' docType methods
#' @rdname profile-methods
setMethods("bandwidth",
  signature = signature(object = "profile"),
  definition = function(object)object@bandwidth
)           

#' readsList
#'
#' @param profile object
#' @return list. A list made off reads objects. One list for each replicate. Must coincide with size of bedfiles
#' @export
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
#' @export
#' docType methods
#' @rdname profile-methods
setMethods("matchList",
  signature = signature(object = "profile"),
  definition = function(object)object@matchList
)
