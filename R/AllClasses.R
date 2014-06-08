
#' profile class
#'
#' Contains all the information necesary for the calculation of profile curves. 
#' @slot regions - List with the regions for which the coverage is going to be calculated
#' @slot bedfiles - Character with the files which contains the reads
#' @slot fragLen - Numeric fragment length to extend the reads
#' @slot bandwidth - Numeric value used to smooth the individual coverages
#' @slot readsList - List of reads objects
#' @slot matchList - List of match objects
#' @exportClass profile
setClass("profile",
  representation(regions = "GRangesList",
                 bedfiles = "character",
                 fragLen = "numeric",
                 bandwidth = "numeric",
                 readsList = "list",
                 matchList = "list"),
  prototype = prototype(regions = GRangesList(),
    bedfiles = "",
    fragLen = 0,
    bandwidth = 1,
    readsList = list(),
    matchList = list()))

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

#' reads class
#'
#' Contains the reads obtained in a ChIP - seq experiment separated by strand (for efficiency) and
#' then by chromosome.
#' @slot reads1 - GRangesList of the reads of the ChIP - Seq experiment that have + strand
#' @slot reads2 - GRangesList of the reads of the ChIP - Seq experiment that have - strand
setClass("reads",
  representation(reads1 = "GRangesList",reads2 = "GRangesList"),
  prototype = prototype(reads1 = GRangesList(),reads2 = GRangesList()))         
