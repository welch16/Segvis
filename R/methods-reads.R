
# Methods for reads class

## Get methods

#' @rdname reads-methods
#' @name reads1
#' @aliases reads
setMethods("reads1",
  signature = signature(object = "reads"),
  definition = function(object)object@reads1
)           

#' @rdname reads-methods
#' @name reads2
#' @aliases reads
setMethods("reads2",
  signature = signature(object = "reads"),
  definition = function(object)object@reads2
)

## Set methods

#' @rdname reads-methods
#' @name setReads1
#' @aliases reads
setMethods("setReads1",
  signature = signature(object = "reads",r1 = "GRangesList"),
  definition = function(object,r1){
    object@reads1 = r1
    return(object)
})    

#' @rdname reads-methods
#' @name setReads2
#' @aliases reads
setMethods("setReads2",
  signature = signature(object = "reads",r2 = "GRangesList"),
  definition = function(object,r2){
    object@reads2 = r2
    return(object)
})

#' @rdname reads-methods
#' @name show
#' @aliases reads
setMethods("show",
  signature = signature(object = "reads"),
  definition = function(object){
    chr = names(reads1(object))
    len1 = sapply(reads1(object),FUN = length)
    len2 = sapply(reads2(object),FUN = length)
    cat("chr","    Reads1","   Reads2 \n")
    cat("---------------------------\n")
    for(c in chr){
      cat(c,ifelse(nchar(c)==4,"  "," "),len1[c]," ",len2[c],"\n")
    }
})

#' @rdname reads-methods
#' @name length
#' @aliases reads
setMethods("length",
  signature = signature(x = "reads"),
  definition = function(x){   
    len1 = sapply(reads1(x),FUN = length)
    len2 = sapply(reads2(x),FUN = length)
    return(sum(len1)+sum(len2))
})           



