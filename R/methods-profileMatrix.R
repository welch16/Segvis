
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
  signature = signature(object = "profileMatrix",newName = "character"),
  definition = function(object,newName){
    object@name <- newName
    return(object)
})    

# @rdname profileMatrix-methods
# @name setRegions
# @aliases profileMatrix
setReplaceMethod("regions",
  signature = signature(object = "profileMatrix",newRegions = "GRanges"),
  definition = function(object,newRegions){
    stopifnot(class(object) == "GRanges")
    object@regions = newRegions
    return(object)    
})

# @rdname profileMatrix-methods
# @name setProfileMat
# @aliases profileMatrix
setReplaceMethod("profileMat",
  signature = signature(object = "profileMatrix",newProfileMat = "matrix"),
  definition = function(object,newProfileMat){
    stopifnot(class(newProfileMat) == "matrix")
    stopifnot(length(regions(object))==nrow(newProfileMat))
    object@profileMat = newProfileMat
    return(object)    
})           

# @rdname profileMatrix-methods
# @name setBandwidth
# @aliases profileMatrix
setReplaceMethod("bandwidth",
  signature = signature(object = "profileMatrix",newBandwidth = "numeric"),
  definition = function(object,newBandwidth){
    stopifnot(newBandwidth >=1)
    stopifnot(newBandwidth %% 2 == 0)
    object@bandwidth = newBandwidth
    return(object)
})           

# @rdname profileMatrix-methods
# @name setNormConst
# @aliases profileMatrix
setReplaceMethod("normConst",
  signature = signature(object = "profileMatrix",newNormConst = "numeric"),
  definition = function(object,newNormConst){
    stopifnot(newNormConst > 0)
    object@normConst = newNormConst
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
# @name 


                   
