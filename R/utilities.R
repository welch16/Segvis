
#' @import data.table
NULL


#' @export
.data.table.GRanges <- function(x)
{
  dt = data.table(seqnames = as.character( seqnames(x)),
    start = start(x),end = end(x),
    strand = as.character(strand(x)))
  return(dt)
}

#' @export
.GRanges.data.table <- function(x)
{
  stopifnot(c("seqnames","start","end","strand") %in% names(x))
  return(GRanges(seqnames = x$seqnames,ranges =
    IRanges(start = x$start,end = x$end),
    strand = x$strand))                 
}

#' @export
.separate.reads <- function(greads,chrom,side,mc)
{
  chr_reads = mclapply(chrom,function(ch,greads,side){
    return(greads[seqnames == ch & strand == side])},
    greads,side,mc.cores = mc)
  return(chr_reads)
}
