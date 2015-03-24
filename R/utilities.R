
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
  return(GRanges(seqnames = x[,(seqnames)],
    ranges = .IRanges.data.table(x),
    strand = x[,(strand)]))             
}

#' @export
.IRanges.data.table <- function(x)
{
  stopifnot(c("seqnames","start","end","strand") %in% names(x))
  return(IRanges(start = x[,(start)],end = x[,(end)]))
}


