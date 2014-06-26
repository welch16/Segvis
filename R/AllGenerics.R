
# Generic methods for several classes

## Get methods

#' name
#'
#' @param profile object
#' @return chracter. The name of the profile object
#' @export
#' docType methods
setGeneric("name",
  function(object,...)
  standardGeneric("name")
)           

#' regions
#'
#' @param profile object
#' @return GRangesList. The regions used to calculate the coverage plots, as a GRangesList separated by chromosome
#' @export
#' @docType methods
setGeneric("regions",
  function(object,...)
  standardGeneric("regions")
)

## Set methods

#' setName
#'
#' @param profile object
#' @param newName character
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setName",
  function(object,newName)
  standardGeneric("setName")           
)

#' setRegions
#'
#' @param profile object
#' @param newRegions GRangesList
#' @return profile object
#' @export
#' docType methods
#' @rdname profile-methods
setGeneric("setRegions",
  function(object,newRegions)
  standardGeneric("setRegions")
)

# Generic methods for profile class

##  Get methods

#' files
#'
#' @param profile object
#' @return character. A chracter vector with the names of the bedfiles used to create the readsList object
#' @export
#' @docType methods
setGeneric("files",
  function(object,...)
  standardGeneric("files")           
)

#' fileFormat
#' @param profile object
#' @return character. A character that contains the file format used
#' @export
#' @docType methods
setGeneric("fileFormat",
  function(object,...)
  standardGeneric("fileFormat")         
)

#' maxBandwidth
#'
#' @param profile object
#' @return numeric. A number with the max bandwidth possible to smooth the profiles
#' @export
#' @docType methods
#' @rdname profile-methods
#' maxBandwidth
#'
setGeneric("maxBandwidth",
  function(object,...)
  standardGeneric("maxBandwidth")           
)
           
#' fragLen
#'
#' @param profile object
#' @return numeric. A numeric value representing the width of the extended fragment reads
#' @export
#' @docType methods
setGeneric("fragLen",
  function(object,...)
  standardGeneric("fragLen")           
)

#' readsList
#'
#' @param profile object
#' @return list. A list made off reads objects. One list for each replicate. Must coincide with size of bedfiles
#' @export
#' @docType methods
setGeneric("readsList",
  function(object,...)
  standardGeneric("readsList")           
)

#' matchList
#'
#' @param profile object
#' @return list. A list made off match objects. One list for each replicate. Must coincide with size of bedfiles
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("matchList",
  function(object,...)
  standardGeneric("matchList")           
)

#' profileCurve
#' @param profile object
#' @return RleList. A list made of an Rle object for each region
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("profileCurve",
  function(object,...)
  standardGeneric("profileCurve")
)  
    
## Set methods

#' setMaxBandwidth
#'
#' @param profile object
#' @param newMaxBandwidth Numeric value, must be odd and greater or equal than one
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setMaxBandwidth",
  function(object,newMaxBandwidth)
  standardGeneric("setMaxBandwidth")
)           

#' setFragLen
#'
#' @param profile object
#' @param newFragLen Numeric value, must be greater of equal to zero
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setFragLen",
  function(object,newFragLen)
  standardGeneric("setFragLen")
)           

#' loadReads
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("loadReads",
  function(object,mc)
  standardGeneric("loadReads")
)           

#' matchReads
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("matchReads",
  function(object,mc)
  standardGeneric("matchReads")
)           

#' getCoverage
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("getCoverage",
  function(object,mc)
  standardGeneric("getCoverage")           
)

#' buildProfileMat
#'
#' @param profile object
#' @param bw, the bandwidth used to smooth the profiles
#' @param mc, the number of cores used with parallel
#' @return list object
#' @export
#' @docType methods
#' rdname methods-profile
setGeneric("buildProfileMat",
  function(object,bw,mc)
  standardGeneric("buildProfileMat")
)

# Generic Methods for reads class

## Get methods

#' reads1
#'
#' @param reads object
#' @return GRangesList The reads of the "+" strand
#' @export
#' @docType methods
#' @rdname reads-methods
setGeneric("reads1",
  function(object,...)
  standardGeneric("reads1")
)           

#' reads2
#'
#' @param reads object
#' @return GRangesList The reads of the "-" strand
#' @export
#' @docType methods
#' @rdname reads-methods
setGeneric("reads2",
  function(object,...)
  standardGeneric("reads2")           
)

## Set methods

#' setReads1
#'
#' @param reads object
#' @param r1 GRangesList object, the new reads to set on the reads object
#' @return reads object
#' @export
#' @docType methods
#' @rdname reads-methods
setGeneric("setReads1",
  function(object,r1)
  standardGeneric("setReads1")
)          

#' setReads2
#'
#' @param reads object
#' @param r2 GRangesList object, the new reads to set on the reads object
#' @return reads object
#' @export
#' @docType methods
#' @rdname reads-methods
setGeneric("setReads2",
  function(object,r2)
  standardGeneric("setReads2")
)          

# Generic Methods for match class

## Get methods

#' match1
#'
#' @param match object
#' @return list The match of the reads with "+" strand
#' @export
#' @docType methods
#' @rdname match-methods
setGeneric("match1",
  function(object,...)
  standardGeneric("match1")
)           

#' match2
#'
#' @param match object
#' @return list The match of the reads with "-" strand
#' @export
#' @docType methods
#' @rdname match-methods
setGeneric("match2",
  function(object,...)
  standardGeneric("match2")
)           

# Generic methods for profileMatrix class

## Get methods

#' profMat
#'
#' @param object. profileMatrix
#' @return matrix object. A matrix with a the profile for each region
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("profMat",
  function(object)
  standardGeneric("profMat")           
)        

#' bandwidth
#'
#' @param object. profileMatrix
#' @return numeric. The bandwidth used to smooth the profile in profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("bandwidth",
  function(object)
  standardGeneric("bandwidth")           
)           

#' normConst
#'
#' @param object. profileMatrix
#' @return numeric. The constant used to normalize the profile
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("normConst",
  function(object)
  standardGeneric("normConst")           
)           

## Set methods

#' setProfMat
#'
#' @param object. profileMatrix
#' @param newProfMat matrix object. A matrix with a the profile for each region
#' @return profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("setProfMat",
  function(object,newProfMatrix)
  standardGeneric("setProfMat")
)           

#' setBandwidth
#'
#' @param object. profileMatrix
#' @param newBandwidth. Numeric value with the bandwidth used to smooth the profile in profileMatrix object
#' @return profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("setBandwidth",
  function(object,newBandwidth)
  standardGeneric("setBandwidth")
)           

#' setNormConst
#'
#' @param object. profileMatrix
#' @param newNormConst Numeric value with the constant used to normalize the profile
#' @return profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("setNormConst",
  function(object,newNormConst)
  standardGeneric("setNormConst")
)           

            
