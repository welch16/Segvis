
#' profile class
#'
#' Contains all the information necessary for the calculation of profile curves.
#' @slot name - Character with the name of the profiles
#' @slot regions - List with the regions for which the coverage is going to be calculated
#' @slot files - Character with the files which contains the reads
#' @slot fileFormat - Character with the file format used for the reads
#' @slot fragLen - Numeric fragment length to extend the reads
#' @slot bandwidth - Numeric value used to smooth the individual coverages
#' @slot readsList - List of reads objects
#' @slot matchList - List of match objects
#' @slot profileCuve - RleList - For each region, there is a Rle object
#' @slot .haveRegions - logical - Indicates if the object have the regions loaded
#' @slot .haveReads - logical - Indicates if the object have the reads loaded
#' @exportClass profile
setClass("profile",
  representation(name = "character",
                 regions = "GRangesList",
                 files = "character",
                 fileFormat = "character",
                 fragLen = "numeric",
                 bandwidth = "numeric",
                 readsList = "list",
                 matchList = "list",
                 profileCurve = "list",
                 .haveRegions = "logical",
                 .haveReads = "logical"
                 ),
  prototype = prototype(name = "",
    regions = GRangesList(),
    files = "",
    fileFormat = "",
    fragLen = 0,
    bandwidth = 1,
    readsList = list(),
    matchList = list(),
    profileCurve = list(),
    .haveRegions = FALSE,
    .haveReads = FALSE)
)    

setValidity("profile",
  function(object){
  # Checks that readsList and matchList have same length
  match.length = length(object@matchList)
  reads.length = length(object@readsList)
  return(match.length == reads.length & object@fragLen >=0 & object@bandwidth >=1 & tolower(object@fileFormat) == "bam")
}
)  

#' match class
#'
#' Contains the match of a set a reads of a ChIP - Seq experiment given a set
#' of regions in the genome. We are obtaining a match for each strand with the
#' match_reads.cpp function
#' @slot match1 - List of positions of the reads in list that match with each region and have + strand
#' @slot match2 - List of positions of the reads in list that match with each region and have - strand
setClass("match",
  representation(match1 = "list",match2 = "list"),
  prototype = prototype(match1 = list(),
    match2 = list()))

setValidity("match",
  function(object) return(length(object@match1) == length(object@match2))
)            

#' reads class
#'
#' Contains the reads obtained in a ChIP - seq experiment separated by strand (for efficiency) and
#' then by chromosome.
#' @slot reads1 - GRangesList of the reads of the ChIP - Seq experiment that have + strand
#' @slot reads2 - GRangesList of the reads of the ChIP - Seq experiment that have - strand
setClass("reads",
  representation(reads1 = "GRangesList",reads2 = "GRangesList"),
  contains = "GRangesList",
  prototype = prototype(reads1 = GRangesList(),reads2 = GRangesList()))         

setValidity("reads",
  function(object)return(length(object@reads1) == length(object@reads2))
)            






