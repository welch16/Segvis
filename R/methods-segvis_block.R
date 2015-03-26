
# Methods for the profileMatrix class

## Get methods

#' @rdname methods-segvis_block-gs
#' @name name
setMethod("name",
  signature = signature(object = "segvis_block"),
  definition = function(object)object@name
)

#' @rdname methods-segvis_block-gs
#' @name regions
setMethod("regions",
  signature = signature(object = "segvis_block"),
  definition = function(object)object@regions
)

#' @rdname methods-segvis_block-gs
#' @name cover_table
setMethod("cover_table",
  signature = signature(object = "segvis_block"),
  definition = function(object)object@cover_table
)           

#' @rdname methods-segvis_block-gs
#' @name bandwidth
setMethod("bandwidth",
  signature = signature(object = "segvis_block"),
  definition = function(object)object@bandwidth
)           

#' @rdname methods-segvis_block-gs
#' @name normConst
setMethod("normConst",
  signature = signature(object = "segvis_block"),
  definition = function(object)object@normConst
)

## Set methods

#' @rdname methods-segvis_block-gs
#' @name name
setReplaceMethod("name",
  signature = signature(object = "segvis_block",value = "character"),
  definition = function(object,value){
    object@name = value
    return(object)
})    

#' @rdname methods-segvis_block-gs
#' @name regions
setReplaceMethod("regions",
  signature = signature(object = "segvis_block",value = "GRanges"),
  definition = function(object,value){    
    stopifnot(class(value) == "GRanges")
    object@regions = value
    return(object)    
})

#' @rdname methods-segvis_block-gs
#' @name cover_table
setReplaceMethod("cover_table",
  signature = signature(object = "segvis_block",value = "data.table"),
  definition = function(object,value){    
    stopifnot(class(value) == "data.table")
    stopifnot(length(regions(object))==nrow(value))
    object@profileMat = value
    return(object)    
})           

#' @rdname methods-segvis_block-gs
#' @name bandwidth
setReplaceMethod("bandwidth",
  signature = signature(object = "segvis_block",value = "numeric"),
  definition = function(object,value){
    stopifnot(value >=1)
    stopifnot(value %% 2 == 1)
    object@bandwidth = value
    return(object)
})           

#' @rdname methods-segvis_block-gs
#' @name normConst
setReplaceMethod("normConst",
  signature = signature(object = "segvis_block",value = "numeric"),
  definition = function(object,value){
    stopifnot(value > 0)
    object@normConst = value
    return(object)
})                        

# @rdname profileMatrix-methods
# @name meanProfile
setMethods("meanProfile",
  signature = signature(object = "profileMatrix",trim = "numeric"),
  definition = function(object,trim){
    if(missing(trim))trim=0
    stopifnot(is.numeric(trim))
    mat = profileMat(object)
    return(sapply(1:ncol(mat),function(i,mat)
      mean(mat[,i],trim = trim,na.rm = TRUE),mat))
})
    
# @rdname profileMatrix-methods
# @name show
setMethod("show",
  signature = signature(object = "segvis_block"),
  definition = function(object){
    cat("Profile matrix for",name(object),"\n")
    cat("Bandwidth:",bandwidth(object),"\n")
    if( object@.isScaled){
      cat("The profile matrix IS scaled\n")
    }else{
      cat("The profile matrix IS NOT scaled\n")
    }
    show(regions(object))
})

# @rdname profileMatrix-methods
# @name normalize.matrix
setMethods("normalize.matrix",
  signature = signature(object = "profileMatrix",value = "numeric",base = "numeric"),
  definition = function(object, value,base){
    if(missing(base))base = 1000000
    if(missing(value))value = base/normConst(object)
    stopifnot(value > 0)
    profileMat(object) = value * profileMat(object)
    object@.isScaled = TRUE
    return(object) 
})
           
# @rdname profileMatrix-methods
# @name subset.pm
setMethod("subset.pm",
  signature = signature(object = "profileMatrix",condition = "ANY"),
  definition = function(object, condition){    
    cond = .subset_profileMat_logical(object,substitute(condition))
    return(.filter_profileMat(object,cond))
})

# @rdname profileMatrix-methods
# @name addColumn
setMethods("addColumn",
  signature = signature(object = "profileMatrix",name = "character",col = "ANY"),
  definition = function(object,name,col){
    stopifnot(length(col) == length(regions(object)))
    elementMetadata(regions(object))@listData[[name]] = col
    return(object)
})


