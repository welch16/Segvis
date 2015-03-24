
#' @import data.table

#' @title .separate.by.chrom
#'
#' @description Separates a set of reads by chromosome while filters by strand (and
#' have the posibility to sort)
#'
#' @param greads A data.table with at least the columns seqnames, start, end and strand
#'
#' @param chrom A vector with all the chromosome used to split the table
#'
#' @param st A character that represents the strand: "+", "-" or "*";
#'  for forward, reverse or any of them, respectively
#'
#' @param mc A numeric value indicating the number of cores
#'
#' @param sort Optional boolean to see if the intervals are going to be sorted
#'
#' @return list with greads separated 
#' 
#' @export
separate.by.chrom <- function(greads,chrom,st,mc,sort=FALSE)
{
  chr_reads = mclapply(chrom,function(ch,greads,st){
    return(greads[seqnames == ch & strand == st])},
    greads,st,mc.cores = mc)
  if(sort){
    if(st == "-"){ # sort by end
      chr_reads = mclapply(chr_reads,function(x)
        return(x[order(x[,(end)])]),mc.cores = mc)                
    }else{ # sort by start
      chr_reads = mclapply(chr_reads,function(x)
        return(x[order(x[,(start)])]),mc.cores = mc)                
    }
  }
  return(chr_reads)
}


.find.overlaps <- function(reads,regions)
{
  ov = findOverlaps(.IRanges.data.table(reads),
    .IRanges.data.table(regions))
  return(ov)
}

.match.reads <- function(reads,overlaps)
{
  reads[,match:= 0L]
  reads[queryHits(overlaps),match:=subjectHits(overlaps)]
  return(reads)
}

