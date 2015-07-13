#' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

#' reads class description
#'
#' Contains the reads obtained in a ChIP - seq experiment separated by strand and then by chromosome. It has one component for each strand which are object of the data.table class with a match column to identify the regions
#'
#' \section{Slots}{
#'   \describe{
#'     \item{\code{readsF}:}{List of data.table objects containing the reads of the ChIP - Seq experiment that have + strand.}
#'     \item{\code{readsR}:}{List of data.table object containning the reads of the ChIP - Seq experiment that have - strand}
#'   }
#' }
#' @seealso \code{\link{loadReads}}
#'
#' @name reads-class
#' @rdname reads-class
#' @exportClass reads
setClass("reads",
  representation(readsF = "list",readsR = "list"),
  prototype = prototype(readsF = list(),readsR = list()))

setValidity("reads",
  function(object){
    return(length(object@readsF) == length(object@readsR))
})            

#' Segvis class description
#'
#' This object is the base class of the segvis package. It contains all the information necessary for the calculation of coverage curves.
#'
#' \section{Slots}{
#'   \describe{
#'     \item{\code{name}:}{Character with the name of the profiles}
#'     \item{\code{regions}:}{GRanges object with the regions for which the profile want to be calcualted}
#'     \item{\code{file}:}{Character with the name of the file that contains the reads}
#'     \item{\code{maxBandwidth}:}{Numeric value with the maximum bandwidth accepted when smoothing profiles. Must be odd}
#'     \item{\code{ fragLen}:}{Numeric value with the fragment length to resize the reads (if it is zero then it doesn't resize the reads)}
#'     \item{\code{chr}:}{character vector, with the chromosomes to be considered}
#'     \item{\code{isPET}:}{logical, Indicates is the reads come from a PET experiments or a SET experiment}
#'     \item{\code{reads}:}{Reads object, which contains the reads used to build the profile separated by strand}
#'     \item{\code{profiles}:}{RleList - For each region, there is a Rle object}
#'     \item{\code{.haveRegions}:}{logical - Indicates if the object have the regions loaded}
#'     \item{\code{.haveReads}:}{logical - Indicates if the object have the reads loaded}
#'     \item{\code{.readsMatched}:}{logical - Indicates if the read have been matched to the regions}
#'     \item{\code{.coverageCalculated}:}{logical - Indicates if the coverage has been calculated]}
#'    }
#' }
#' @seealso \code{\link{buildSegvis}}
#' @name segvis-class
#' @rdname segvis-class
#' @exportClass segvis
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

setValidity("segvis",
  function(object){
  return(object@fragLen >=0 & object@maxBandwidth >=1)
})  

#' segvis_block class description
#'
#' Contains a data.table with the individual coverages for each region
#'
#' \section{Slots}{
#'   \describe{
#'    \item{\code{name}:}{Character with the name of the profiles}
#'    \item{\code{regions}:}{GRanges object with the regions of the individuals profiles}
#'    \item{\code{cover}:}{table data.table with the individuals coverages for each region}
#'    \item{\code{bandwidth}:}{Numeric value used to smooth the individual profiles}
#'    \item{\code{normCont}:}{Numeric normalizing contanst}
#'    \item{\code{.isScaled}:}{Logical value representing if the profiles are scaled to \code{normConst}}
#'    }
#' }
#' @seealso \code{\link{Segvis_block}}
#' @exportClass segvis_block
#' @name segvis_block-class
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

#' segvis_block_list class description
#'
#' Contains a list of \code{segvis_block} objects
#'
#' \section{Slots}{
#'   \describe{
#'    \item{\code{list}:}{A list of \code{segvis_block} objects or several \code{segvis_block} objects separated by commas}
#'    }
#' }
#' @exportClass segvis_block_list
#' @name segvis_block_list-class
#' @rdname segvis_block_list-class
setClass("segvis_block_list",         
  prototype = prototype(elementType = "segvis_block"),
  contains = "list"
)      


