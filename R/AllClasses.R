#' @title reads class description
#'
#' @description Contains the reads obtained in a ChIP - seq experiment
# separated by strand and then by chromosome.
#'
#' @slot readsF List of data.table objects containing the reads of the
#' ChIP - Seq experiment that have + strand.
#'
#' @slot readsR List of data.table object containning the reads of the
#' ChIP - Seq experiment that have - strand
#'
#' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{loadReads}}
setClass("reads",
  representation(readsF = "list",readsR = "list"),
  prototype = prototype(readsF = list(),readsR = list()))

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
#' @slot chr character vector, with the chromosomes to be considered
#'
#' @slot isPET logical, Indicates is the reads come from a PET experiments or a SET experiment
#'
#' @slot reads Reads object, which contains the reads used to build the profile separated by strand
#'
#' @slot profiles RleList - For each region, there is a Rle object
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
#'
#' @rdname segvis-class
setClass("segvis",
  representation(name = "character",
                 regions = "GRanges",
                 file = "character",
                 maxBandwidth = "numeric",
                 fragLen = "numeric",
                 chr = "character",
                 isPET = "logical",                 
                 reads = "reads",
                 profiles = "list",
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
                 chr = "",
                 isPET = FALSE,
                 reads = new("reads"),
                 profiles = list(),
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
  return(object@fragLen >=0 & object@maxBandwidth >=1)
})  

#' @title segvis_block class description
#'
#' @description Contains a data.table with the individual coverages for each region
#'
#' @slot name Character with the name of the profiles
#'
#' @slot regions GRanges object with the regions of the individuals profiles
#'
#' @slot cover_table data.table with the individuals coverages for each region
#'
#' @slot bandwidth Numeric value used to smooth the individual profiles
#'
#' @slot normCont Numeric normalizing contanst
#'
#' @slot .isScaled Logical value representing if the profiles are scaled to \code{normConst}
#'
#' @exportClass segvis_block
#'
#' @rdname segvis_block-class
setClass("segvis_block",
  representation(name = "character",
                 regions = "GRanges",
                 cover_table = "data.table",
                 bandwidth = "numeric",
                 normConst = "numeric",
                 .isScaled = "logical"),
  prototype = prototype(name = "",
                 regions = GRanges(),
                 cover_table = data.table(),
                 bandwidth = 1,
                 normConst = 1,
                 .isScaled = FALSE)
)

#' @title segvis_block_list class description
#' @description Contains a list of segvis objects
#' @exportClass segvis_block_list
setClass("segvis_block_list",         
  prototype = prototype(elementType = "segvis_block"),
  contains = "list"
)      


