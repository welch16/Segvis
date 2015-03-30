
# Methods for reads class

## Get methods

#' @rdname methods-reads-gs
#' @name readsF
setMethod("readsF",
  signature = signature(object = "reads"),
  definition = function(object)object@readsF
)

#' @rdname methods-reads-gs
#' @name readsR
setMethod("readsR",
  signature = signature(object = "reads"),
  definition = function(object)object@readsR
)

## Set methods

#' @rdname methods-reads-gs
#' @name readsF
setReplaceMethod("readsF",
  signature = signature(object = "reads",value = "list"),
  definition = function(object,value){
    object@readsF = value
    return(object)
})

#' @rdname methods-reads-gs
#' @name readsR
setReplaceMethod("readsR",
  signature = signature(object = "reads",value = "list"),
  definition = function(object,value){
    object@readsR = value
    return(object)
})
