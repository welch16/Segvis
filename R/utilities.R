
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

localMinima <- function(x)
{
  # Use -Inf instead if x is numeric (non-integer)
  y = diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y = cumsum(rle(y)$lengths)
  y = y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}

localMaxima <- function(x) localMinima(-x)
