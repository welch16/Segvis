##' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

# .ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)

##' @rdname SegvisData
##' @export
setClass("SegvisData",
         contains = "GRanges",
         representation = representation(
           files = "character",
           is_pet = "logical",
           frag_len = "numeric",
           covers = "list",
           fwd_covers = "list",
           bwd_covers = "list",
           nreads = "numeric"
           ),
         prototype = prototype(
           files = "",
           is_pet = FALSE,
           frag_len = 1L,
           covers = list(),
           fwd_covers = list(),
           bwd_covers = list(),
           nreads = 0L
           ))

setValidity("SegvisData",
            function(object){
              return(length(object@files) >= 1 &
                       all(object@frag_len >= 0))
            })

##' SegvisData object and constructors
##'
##' \code{SegvisData} is a subclass of \code{GenomicRanges}, used to visualize
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
##' @param mc.cores a numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##'
##' @return  \code{SegvisData} returns a \code{SegvisData} object which contains
##' the forward, backward and aggregated coverage for every aligned read file in
##' \code{files}.
##' @aliases SegvisData SegvisData-class
##'
##' @docType class
##
## @examples
##
## dr = system.file("extdata","example",package = "Segvis",mustWork = TRUE)
## files = list.files(dr,pattern = "bam",full.names =TRUE)
## files = files[grep("bai",files,invert = TRUE)]
## reg = list.files(dr,pattern = "narrow",full.names =TRUE)
## reg = readBedFile(reg[1])
## segvis = SegvisData(regions = reg,files)
##
##' @rdname SegvisData
##' @export
SegvisData = function(regions,files,is_pet = rep(FALSE,length(files)),
                       frag_len = 1,mc.cores = getOption("mc.cores" , 2L))
{
  stopifnot(is.character(files),all(file.exists(files)),
            is.logical(is_pet),is.numeric(frag_len))
  if(length(is_pet) == 1){
    is_pet = rep(is_pet,length(files))
  }
  if(length(frag_len) == 1){
    frag_len = rep(frag_len,length(files))
  }
  stopifnot(length(files) == length(is_pet) ,
            length(files) == length(frag_len))
  
  if(Sys.info()[["sysname"]] == "Windows"){
    parallel_param = SnowParam(workes = mc.cores , type = "SOCK")
  }else{
    parallel_param = MulticoreParam(workers = mc.cores)
  }
  
  fwd_covers = bpmapply(.readFileCover,
                        files,
                        is_pet,
                        frag_len,
                        MoreArgs = list(st = "+"),
                        BPPARAM = parallel_param,
                        SIMPLIFY = FALSE)
  
  bwd_covers = bpmapply(.readFileCover,
                        files,
                        is_pet,
                        frag_len,
                        MoreArgs = list(st = "-"),
                        BPPARAM = parallel_param,
                        SIMPLIFY = FALSE)
  
  covers = mapply("+",fwd_covers,bwd_covers,
                  SIMPLIFY = FALSE)
  nreads = bpmapply(.countReads,
                    files,is_pet,
                    BPPARAM = parallel_param,
                    SIMPLIFY = TRUE)
  names(fwd_covers) = NULL
  names(bwd_covers) = NULL
  names(covers) = NULL
  names(nreads) = NULL
  frag_len[is_pet] = 0
  new("SegvisData",regions,
      files = files,
      is_pet = is_pet,
      frag_len = frag_len,
      covers = covers,
      fwd_covers = fwd_covers,
      bwd_covers = bwd_covers,
      nreads = nreads
  )
}
