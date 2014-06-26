
# Methods for the profileMatrix class

## Get methods

#' name
#'
#' @param object profileMatrix
#' @return character. The name of the object
#' @docType methods
#' @rdname profileMatrix-methods
setMethods("name",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@name
)

#' regions
#'
#' @param object. profileMatrix
#' @return GRanges object. The list of regions for which the profileMatrix was calculated
#' @docType methods
#' @rdname profileMatrix-methods
setMethods("regions",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@regions
)

#' profMat
#'
#' @param object. profileMatrix
#' @return matrix object. A matrix with a the profile for each region
#' @docType methods
#' @rdname profileMatrix-methods
setMethods("profMat",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@profMat
)           

#' bandwidth
#'
#' @param object. profileMatrix
#' @return numeric. The bandwidth used to smooth the profile in profileMatrix object
#' @docType methods
#' @rdname profileMatrix-methods
setMethods("bandwidth",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@bandwidth
)           

#' normConst
#'
#' @param object. profileMatrix
#' @return numeric. The constant used to normalize the profile
#' @docType methods
#' @rdname profileMatrix-methods
setMethods("normConst",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@normConst
)           





              

