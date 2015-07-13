
# Methods for reads class

#### readsF

##' @name readsF
##' @rdname readsF
##' @docType methods
##' @aliases readsF
##' @exportMethod readsF
setMethod("readsF",
    signature = signature(object = "reads"),
    definition = function(object)object@readsF
)

##' @name readsF
##' @rdname readsF
##' @docType methods
##' @aliases readsF<-
##' @exportMethod readsF<-
setReplaceMethod("readsF",
    signature = signature(object = "reads",value = "list"),
    definition = function(object,value){
        object@readsF <- value
        return(object)
})

#### readsR

##' @name readsR
##' @rdname readsR
##' @docType methods
##' @aliases readsR
##' @exportMethod readsR
setMethod("readsR",
    signature = signature(object = "reads"),
    definition = function(object)object@readsR
)


##' @name readsR
##' @rdname readsR
##' @docType methods
##' @aliases readsR<-
##' @exportMethod readsR<-
setReplaceMethod("readsR",
    signature = signature(object = "reads",value = "list"),
    definition = function(object,value){
        object@readsR <- value
        return(object)
})
