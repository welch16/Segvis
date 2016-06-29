##' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

##' @rdname SegvizData
##' @export
setClass("SegvizData",
         contains = "GRanges",
         representation = representation(
           files = "character",
           is_pet = "logical",
           frag_len = "numeric",
           covers = "list",
           fwd_covers = "list",
           bwd_covers = "list"
           ),
         prototype = prototype(
           files = "",
           is_pet = FALSE,
           frag_len = 1L,
           covers = list(),
           fwd_covers = list(),
           bwd_covers = list()
           ))

setValidity("SegvizData",
            function(object){
              return(length(object@files) >= 1 &
                       all(object@frag_len >= 1))
            })

##' SegvizData object and constructors
##'
##' \code{SegvizData} is a subclass of \code{GenomicRanges}, used to visualize
##' high-thoughput sequencing experiments across a set of user - defined
##' genomic regions.
##'
##' @param regions a \code{GRanges} object with the regions to be considered.
##' @param files a character vector with the location of the bam files that
##' contain the aligned reads.
##' @param is_pet a logical vector with the same length as files indicating if
##' the aligned reads in the respective bam file are paired.
##' @param frag_len a numeric vector representing the average fragment length
##' to which the aligned reads in their respective bam file are going to be
##' extended. For PE reads, this parameter is not considered.
##' @param mc.cores A numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##'
##' @return A SegvizData object.
##' @aliases SegvizData SegvizData-class
##'
##' @docType class
##'
##' @examples
##' a = 1
##' @rdname SegvizData
##' @export
SegvizData <- function(regions,files,is_pet = rep(FALSE,length(files)),
                       frag_len = 1,mc.cores = getOption("mc.cores" , 2L))
{
  stopifnot(is.character(files),is.logical(is_pet),is.numeric(frag_len))
  if(length(is_pet) == 1){
    is_pet = rep(is_pet,length(files))
  }
  if(length(frag_len) == 1){
    frag_len = rep(frag_len,length(files))
  }
  stopifnot(length(files) == length(is_pet) ,
            length(files) == length(frag_len))
  fwd_covers = mcmapply(.readFileCover,
                        files,
                        is_pet,
                        frag_len,
                        MoreArgs = list(st = "+"),
                        mc.cores = mc.cores,
                        SIMPLIFY = FALSE)
  bwd_covers = mcmapply(.readFileCover,
                        files,
                        is_pet,
                        frag_len,
                        MoreArgs = list(st = "-"),
                        mc.cores = mc.cores,
                        SIMPLIFY = FALSE)
  covers = mapply("+",fwd_covers,bwd_covers,
                  SIMPLIFY = FALSE)
  new("SegvizData",regions,
      files = files,
      is_pet = is_pet,
      frag_len = frag_len,
      covers = covers,
      fwd_covers = fwd_covers,
      bwd_covers = bwd_covers
  )
}






##' reads class description
##'
##' Contains the reads obtained in a ChIP - seq experiment separated by strand and then by chromosome. It has one component for each strand which are object of the data.table class with a match column to identify the regions
##'
##' @slot readsF List of data.table objects containing the reads of the ChIP - Seq experiment that have + strand.
##' @slot readsR List of data.table object containning the reads of the ChIP - Seq experiment that have - strand
##' @seealso \code{\link{loadReads}}
##'
##' @name reads-class
##' @rdname reads-class
##' @exportClass reads
setClass("reads",
  representation(readsF = "list",readsR = "list"),
  prototype = prototype(readsF = list(),readsR = list()))

setValidity("reads",
  function(object){
    return(length(object@readsF) == length(object@readsR))
})

##' Segvis class description
##'
##' This object is the base class of the segvis package. It contains all the information necessary for the calculation of coverage curves.
##'
##' @slot name Character with the name of the profiles
##' @slot regions GRanges object with the regions for which the profile want to be calcualted
#' @slot file Character with the name of the file that contains the reads
##' @slot maxBandwidth Numeric value with the maximum bandwidth accepted when smoothing profiles. Must be odd
##' @slot fragLen Numeric value with the fragment length to resize the reads (if it is zero then it doesn't resize the reads)
##' @slot chr character vector, with the chromosomes to be considered
##' @slot isPET logical, Indicates is the reads come from a PET experiments or a SET experiment
##' @slot reads Reads object, which contains the reads used to build the profile separated by strand
##' @slot profiles RleList - For each region, there is a Rle object
##' @slot .haveRegions logical - Indicates if the object have the regions loaded
##' @slot .haveReads logical - Indicates if the object have the reads loaded
##' @slot .readsMatched logical - Indicates if the read have been matched to the regions
##' @slot .coverageCalculated logical - Indicates if the coverage has been calculated]
##' @seealso \code{\link{buildSegvis}}
##' @name segvis-class
##' @rdname segvis-class
##' @exportClass segvis
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

##' segvis_block class description
##'
##' Contains a data.table with the individual coverages for each region
##'
##' @slot name Character with the name of the profiles
##' @slot regions GRanges object with the regions of the individuals profiles
##' @slot cover_table data.table with the individuals coverages for each region
##' @slot bandwidth Numeric value used to smooth the individual profiles
##' @slot normConst Numeric normalizing contanst
##' @slot .isScaled Logical value representing if the profiles are scaled to \code{normConst}
##' @seealso \code{\link{Segvis_block}}
##' @exportClass segvis_block
##' @name segvis_block-class
##' @rdname segvis_block-class
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

##' segvis_block_list class description
##'
##' Contains a list of \code{segvis_block} objects
##'
##' @exportClass segvis_block_list
##' @name segvis_block_list-class
##' @rdname segvis_block_list-class
setClass("segvis_block_list",
  prototype = prototype(elementType = "segvis_block"),
  contains = "list"
)


