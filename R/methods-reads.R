
# Methods for reads class

## Get methods

# @rdname reads-methods
# @name readsF
# @aliases reads
setMethod("readsF",
  signature = signature(object = "reads"),
  definition = function(object)object@readsF
)           

# @rdname reads-methods
# @name readsR
# @aliases reads
setMethod("readsR",
  signature = signature(object = "reads"),
  definition = function(object)object@readsR
)

## ## Set methods

## # @rdname reads-methods
## # @name setReads1
## # @aliases reads
## setReplaceMethod("readsF",
##   signature = signature(object = "reads",value = "GRangesList"),
##   definition = function(object,value){
##     object@readsF = value
##     return(object)
## })    

## # @rdname reads-methods
## # @name setReads2
## # @aliases reads
## setReplaceMethod("readsR",
##   signature = signature(object = "reads",value = "GRangesList"),
##   definition = function(object,value){
##     object@readsR = value
##     return(object)
## })

## # @rdname reads-methods
## # @name show
## # @aliases reads
## setMethods("show",
##   signature = signature(object = "reads"),
##   definition = function(object){
##     chr = names(reads1(object))
##     len1 = sapply(readsF(object),FUN = length)
##     len2 = sapply(readsR(object),FUN = length)
##     cat("chr","    Reads1","   Reads2 \n")
##     cat("---------------------------\n")
##     for(c in chr){
##       cat(c,ifelse(nchar(c)==4,"  "," "),len1[c]," ",len2[c],"\n")
##     }
## })




