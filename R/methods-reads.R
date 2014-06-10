
# Methods for reads class

## Get methods

#' reads1
#' @param reads object
#' @return GRangesList The reads of the "+" strand
#' @docType methods
#' @rdname reads-methods
setMethods("reads1",
  signature = signature(object = "reads"),
  definition = function(object)object@reads1
)           

#' reads2
#' @param reads object
#' @return GRangesList The reads of the "-" strand
#' @docType methods
#' @rdname reads-methods
setMethods("reads2",
  signature = signature(object = "reads"),
  definition = function(object)object@reads2
)

## Set methods

#' setReads1
#' @param reads object
#' @param r1 GRangesList object, the new reads to set on the reads object
#' @return reads object
#' @docType methods
#' @rdname reads-methods
setMethods("setReads1",
  signature = signature(object = "reads",r1 = "GRangesList"),
  definition = function(object,r1){
    object@reads1 = r1
    return(object)
})    

#' setReads2
#' @param reads object
#' @param r2 GRangesList object, the new reads to set on the reads object
#' @return reads object
#' @docType methods
#' @rdname reads-methods
setMethods("setReads2",
  signature = signature(object = "reads",r2 = "GRangesList"),
  definition = function(object,r2){
    object@reads2 = r2
    return(object)
})







