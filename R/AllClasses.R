
#' @title match class description
#' @description Contains the match of a set a reads of a ChIP - Seq experiment given a set of regions in the genome. 
#' @slot match1 - List of positions of the reads in list that match with each region and have + strand
#' @slot match2 - List of positions of the reads in list that match with each region and have - strand
#' @seealso \code{\link{match1}} and  \code{\link{match2}}
setClass("match",
  representation(match1 = "list",match2 = "list"),
  contains = "list",         
  prototype = prototype(match1 = list(),
    match2 = list()))

setValidity("match",
  function(object) return(length(object@match1) == length(object@match2))
)            

#' @title reads class description
#' @description Contains the reads obtained in a ChIP - seq experiment separated by strand and then by chromosome.
#' @slot reads1 - GRangesList of the reads of the ChIP - Seq experiment that have + strand
#' @slot reads2 - GRangesList of the reads of the ChIP - Seq experiment that have - strand
#' @seealso \code{\link{reads1}}, \code{\link{reads2}} and loadReads
setClass("reads",
  representation(reads1 = "GRangesList",reads2 = "GRangesList"),
  contains = "GRangesList",
  prototype = prototype(reads1 = GRangesList(),reads2 = GRangesList()))

setValidity("reads",
  function(object)return(length(object@reads1) == length(object@reads2))
)            

#' @title profile class description
#' @description Contains all the information necessary for the calculation of profile curves.
#' @slot name - Character with the name of the profiles
#' @slot regions - GRanges object with the regions for which the profile want to be calcualted
#' @slot file - Character with the name of the file that contains the reads
#' @slot fileFormat - Character with the file format used for the reads
#' @slot maxBandwidth - Numeric - maximum bandwidth accepted when smoothing profiles. Must be odd
#' @slot fragLen - Numeric fragment length to extend the reads
#' @slot remChr - Character - Vector with the chromosomes to be ignored
#' @slot reads - Reads object, which contains the reads used to build the profile separated by strand
#' @slot match - Match object, which contains the 
#' @slot profileCuve - RleList - For each region, there is a Rle object
#' @slot .haveRegions - logical - Indicates if the object have the regions loaded
#' @slot .haveReads - logical - Indicates if the object have the reads loaded
#' @slot .readsMatched - logical - Indicates if the read have been matched to the regions
#' @slot .coverageCalculated - logical - Indicates if the coverage has been calculated]
#' @seealso \code{\link{Profile}}
#' @exportClass profile
setClass("profile",
  representation(name = "character",
                 regions = "GRanges",
                 file = "character",
                 fileFormat = "character",
                 maxBandwidth = "numeric",
                 fragLen = "numeric",
                 remChr = "character",
                 reads = "reads",
                 match = "match",
                 profileCurve = "list",
                 .haveRegions = "logical",
                 .haveReads = "logical",
                 .readsMatched = "logical",
                 .coverageCalculated = "logical"
                 ),
  contains = c("reads","match"),
  prototype = prototype(name = "",
    regions = GRanges(),
    file = "",
    fileFormat = "",
    maxBandwidth = 1,
    fragLen = 0,
    remChr = "",
    reads = new("reads"),
    match = new("match"),
    profileCurve = list(),
    .haveRegions = FALSE,
    .haveReads = FALSE,
    .readsMatched = FALSE,
    .coverageCalculated = FALSE)
)    

setValidity("profile",
  function(object){
  # Checks that readsList and matchList have same length  
  return(object@fragLen >=0&object@maxBandwidth >=1  & tolower(object@fileFormat) == "bam")
}
)  

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
  contains = "GRanges",
  prototype = prototype(name = "",
                 regions = GRanges(),
                 profileMat = matrix(nrow=0,ncol = 0),
                 bandwidth = 1,
                 normConst = 1,
                 .isScaled = FALSE)
)    

