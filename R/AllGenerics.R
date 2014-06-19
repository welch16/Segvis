
# Generic methods for profile class

##  Get methods

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

#' bandwidth
#'
#' @param profile object
#' @return numeric. A numeric value representing the bandwidth used to smooth the average coverage plot
#' @export
#' @docType methods
setGeneric("bandwidth",
  function(object,...)
  standardGeneric("bandwidth")           
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

#' setBandwidth
#'
#' @param profile object
#' @param newBandwidth Numeric value, must be greater of equal to one
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setBandwidth",
  function(object,newBandwidth)
  standardGeneric("setBandwidth")
)           

#' loadReads
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' rdname profile-methods
setGeneric("loadReads",
  function(object,mc)
  standardGeneric("loadReads")
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
