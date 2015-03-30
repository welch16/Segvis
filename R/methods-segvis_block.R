
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

#' @rdname methods-segvis_block-summarize
#' @name summarize
setMethods("summarize",
  signature = signature(object = "segvis_block",FUN = "function",... = "ANY"),
  definition = function(object,FUN,...){
    ## check length of matches
    lengths = cover_table(object)[,length(coord),by = .(chr,match)]
    if(length(u <- unique(lengths[,(V1)])) > 1){
      stop("All regions must have the same length")
    }
    out = copy(cover_table(object))
    out[,center:=0L]
    out[,center:= out[,coord - min(coord) + 1,by = .(chr,match)][,(V1)]]           
    summary = out[,FUN(tagCounts,...),by =.(chr,center)]
    return(summary[,(V1)])
})
                     
# @rdname profileMatrix-methods
# @name show
setMethod("show",
  signature = signature(object = "segvis_block"),
  definition = function(object){
    cat("Segvis profile:",name(object),"\n")
    cat("Bandwidth:",bandwidth(object),"\n")
    if( object@.isScaled){
      cat("The profile matrix IS scaled\n")
    }else{
      cat("The profile matrix IS NOT scaled\n")
    }
    show(regions(object))
})

#' @rdname methods-segvis_block-normalize
#' @name normalize
setMethods("normalize",
  signature = signature(object = "segvis_block",value = "numeric",base = "numeric"),
  definition = function(object, value,base){
    if(missing(base))base = 1000000
    if(missing(value))value = base/normConst(object)
    stopifnot(value > 0)
    cover_table(object)[,tagCounts := value * tagCounts]
    object@.isScaled = TRUE
    return(object) 
})
           
#' @rdname methods-segvis_block-subset
#' @name subset
setMethod("subset",
  signature = signature(object = "segvis_block",condition = "ANY"),
  definition = function(object, condition){
    cond = .subset_logical(object,substitute(condition))
    return(.filter_sb(object,cond))
})

#' @rdname methods-segvis_block-addColumn
#' @name addColumn
setMethods("addColumn",
  signature = signature(object = "segvis_block",name = "character",col = "ANY"),
  definition = function(object,name,col){
    stopifnot(length(col) == length(regions(object)))    
    elementMetadata(regions(object))@listData[[name]] = col
    return(object)
})


