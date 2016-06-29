

##' @rdname files-methods
##' @aliases files
##' @docType methods
##' @exportMethod files
setMethod("files",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@files)

##' @rdname is_pet-methods
##' @aliases is_pet
##' @docType methods
##' @exportMethod is_pet
setMethod("is_pet",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@is_pet)

##' @rdname frag_len-methods
##' @aliases frag_len
##' @docType methods
##' @exportMethod frag_len
setMethod("frag_len",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@frag_len)

##' @rdname covers-methods
##' @aliases covers
##' @docType methods
##' @exportMethod covers
setMethod("covers",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@covers)

##' @rdname fwd_covers-methods
##' @aliases fwd_covers
##' @docType methods
##' @exportMethod fwd_covers
setMethod("fwd_covers",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@fwd_covers)

##' @rdname bwd_covers-methods
##' @aliases bwd_covers
##' @docType methods
##' @exportMethod bwd_covers
setMethod("bwd_covers",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@bwd_covers)

##' @rdname nreads-methods
##' @aliases nreads
##' @docType methods
##' @exportMethod nreads
setMethod("nreads",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@nreads)

