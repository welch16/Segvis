
# Methods for the profileMatrix class

## Get methods

# @rdname profileMatrix-methods
# @name name
# @aliases profileMatrix
setMethod("name",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@name
)

# @rdname profileMatrix-methods
# @name regions
# @aliases profileMatrix
setMethod("regions",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@regions
)

# @rdname profileMatrix-methods
# @name profMat
# @aliases profileMatrix
setMethod("profileMat",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@profileMat
)           

# @rdname profileMatrix-methods
# @name bandwidth
# @aliases profileMatrix
setMethod("bandwidth",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@bandwidth
)           

# @rdname profileMatrix-methods
# @name normConst
# @aliases profileMatrix
setMethod("normConst",
  signature = signature(object = "profileMatrix"),
  definition = function(object)object@normConst
)

## Set methods

# @rdname profileMatrix-methods
# @name name
# @aliases profileMatrix
setReplaceMethod("name",
  signature = signature(object = "profileMatrix",value = "character"),
  definition = function(object,value){
    object@name = value
    return(object)
})    

# @rdname profileMatrix-methods
# @name setRegions
# @aliases profileMatrix
setReplaceMethod("regions",
  signature = signature(object = "profileMatrix",value = "GRanges"),
  definition = function(object,value){    
    stopifnot(class(value) == "GRanges")
    object@regions = value
    return(object)    
})

# @rdname profileMatrix-methods
# @name setProfileMat
# @aliases profileMatrix
setReplaceMethod("profileMat",
  signature = signature(object = "profileMatrix",value = "matrix"),
  definition = function(object,value){    
    stopifnot(class(value) == "matrix")
    stopifnot(length(regions(object))==nrow(value))
    object@profileMat = value
    return(object)    
})           

# @rdname profileMatrix-methods
# @name setBandwidth
# @aliases profileMatrix
setReplaceMethod("bandwidth",
  signature = signature(object = "profileMatrix",value = "numeric"),
  definition = function(object,value){
    stopifnot(value >=1)
    stopifnot(value %% 2 == 1)
    object@bandwidth = value
    return(object)
})           

# @rdname profileMatrix-methods
# @name setNormConst
# @aliases profileMatrix
setReplaceMethod("normConst",
  signature = signature(object = "profileMatrix",value = "numeric"),
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
  signature = signature(object = "profileMatrix"),
  definition = function(object){
    cat("Profile matrix for",name(object),"\n")
    cat("Bandwidth:",bandwidth(object),"\n")
    show(regions(object))
})

# @rdname profileMatrix-methods
# @name subset.pm
setMethod("subset.pm",
  signature = signature(object = "profileMatrix",subset = "ANY"),
  definition = function(object, subset){
    env = list2env(as(regions(object),"data.frame"),parent = parent.frame())
    cond = as.logical(eval(substitute(subset),env))   
    return(new("profileMatrix",name = name(object),regions = regions(object)[cond],
      bandwidth = bandwidth(object),normConst = normConst(object),profileMat = profileMat(object)[cond,],
               .isScaled = object@.isScaled))
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
