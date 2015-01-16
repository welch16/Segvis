
#' @title match class description
#' @description Contains the match of a set a reads of a ChIP - Seq experiment given a set of regions in the genome. 
#' @slot matchF - List of positions of the reads in list that match with each region and have + strand
#' @slot matchR - List of positions of the reads in list that match with each region and have - strand
#' @seealso \code{\link{matchF}} and  \code{\link{matchR}}
setClass("match",
  representation(matchF = "list",matchR = "list"),
  contains = "list",         
  prototype = prototype(matchF = list(),matchR = list()))

setValidity("match",
  function(object){
    return(length(object@matchF) == length(object@matchR))
})            

#' @title reads class description
#' @description Contains the reads obtained in a ChIP - seq experiment separated by strand and then by chromosome.
#' @slot readsF - GRangesList of the reads of the ChIP - Seq experiment that have + strand
#' @slot readsR - GRangesList of the reads of the ChIP - Seq experiment that have - strand
#' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{loadReads}}
setClass("reads",
  representation(readsF = "GRangesList",readsR = "GRangesList"),
  contains = "GRangesList",
  prototype = prototype(readsF = GRangesList(),readsR = GRangesList()))

setValidity("reads",
  function(object){
    return(length(object@readsF) == length(object@readsR))
})            

#' @title segvis class description
#'
#' @description This object is the base class of the segvis package. It contains all the information necessary for the calculation of coverage curves
#' 
#' @slot name Character with the name of the profiles
#'
#' @slot regions GRanges object with the regions for which the profile want to be calcualted
#'
#' @slot file Character with the name of the file that contains the reads
#'
#' @slot maxBandwidth Numeric value with the maximum bandwidth accepted when smoothing profiles. Must be odd
#' 
#' @slot fragLen Numeric value with the fragment length to resize the reads (if it is zero then it doesn't resize the reads)
#'  
#' @slot reads Reads object, which contains the reads used to build the profile separated by strand
#'
#' @slot match Match object, which contains the
#'
#' @slot profileCuve RleList - For each region, there is a Rle object
#'
#' @slot .haveRegions logical - Indicates if the object have the regions loaded
#'
#' @slot .haveReads logical - Indicates if the object have the reads loaded
#'
#' @slot .readsMatched logical - Indicates if the read have been matched to the regions
#'
#' @slot .coverageCalculated logical - Indicates if the coverage has been calculated]
#'
#' @seealso \code{\link{Segvis}}
#'
#' @exportClass segvis
setClass("segvis",
  representation(name = "character",
                 regions = "GRanges",
                 file = "character",
                 maxBandwidth = "numeric",
                 fragLen = "numeric",
                 reads = "reads",
                 match = "match",
                 profileCurve = "list",
                 .haveRegions = "logical",
                 .haveReads = "logical",
                 .readsMatched = "logical",
                 .coverageCalculated = "logical"
                 ),
  prototype = prototype(name = "",
                 regions = GRanges(),
                 file = "",
                 maxBandwidth = 1,
                 fragLen = 0,     
                 reads = new("reads"),
                 match = new("match"),
                 profileCurve = list(),
                 .haveRegions = FALSE,
                 .haveReads = FALSE,
                 .readsMatched = FALSE,
                 .coverageCalculated = FALSE)
)
     # Slots removed:
     # @slot fileFormat - Character with the file format used for the reads
     # @slot remChr - Character - Vector with the chromosomes to be ignored


setValidity("segvis",
  function(object){
  # Checks that readsList and matchList have same length  
  return(object@fragLen >=0 & object@maxBandwidth >=1)
})  

#' @title profileMatrix class description
#' @description Contains a matrix with the individual coverage for each region
#' @slot name - Character with the name of the profiles
#' @slot regions - GRanges with the regions
#' @slot profileMat - Matrix - Actual profile matrix
#' @slot bandwidth - Numeric value used to smooth the individual coverages
#' @slot normConst - Numeric Normalizing constant
#' @slot .isScaled - Logical representing if the profile matrix is scaled
#' @exportClass profileMatrix
setClass("profileMatrix",
  representation(name = "character",
                 regions = "GRanges",
                 profileMat = "matrix",
                 bandwidth = "numeric",                 
                 normConst = "numeric",
                 .isScaled = "logical"),
#  contains = "GRanges",
  prototype = prototype(name = "",
                 regions = GRanges(),
                 profileMat = matrix(nrow=0,ncol = 0),
                 bandwidth = 1,
                 normConst = 1,
                 .isScaled = FALSE)
)    

#' @title profileMatrixList class description
#' @description Contains a list of profileMatrix objects
#' @exportClass profileMatrixList
setClass("profileMatrixList",
  prototype = prototype(elementType = "profileMatrix"),
  contains = "list"
)         

