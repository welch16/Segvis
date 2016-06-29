
##' @importFrom Rsamtools scanBamFlag ScanBamParam countBam
##' @importFrom GenomicAlignments readGAlignments
NULL

.readFileCover <- function(file,is_pet,frag_len,st = "*")
{
  stopifnot(st %in% c("+","-","*"))
  if(is_pet){
    pet_flag = scanBamFlag(isPaired = TRUE,isMinusStrand = (st == "-"))
    param = ScanBamParam(flag = pet_flag)
    cover = coverage(file,param  = param)
  }else{
    if(st != "*"){
      flag = scanBamFlag(isMinusStrand = (st == "-"))
      par = ScanBamParam(flag = flag)
    }else{
      par = NULL
    }
    reads = readGAlignments(file,param = par)
    reads = resize(as(reads,"GRanges"),frag_len)
    cover = coverage(reads)
  }
  return(cover)
}

.countReads <- function(file,is_pet)
{
  if(is_pet){
    pet_flag = scanBamFlag(isPaired = TRUE,isMinusStrand = (st == "-"))
    par = ScanBamParam(flag = pet_flag)
  }else{
    par = ScanBamParam()
  }
  nr = countBam(file,param = par)$records
  if(is_pet){
    nr = nr/2
  }
  return(nr)
}

.dt2gr <- function(x){
  if(any(!names(x)[1:3] %in% c("seqnames","start","end"))){
    x = x[,1:3,with = FALSE]
    setnames(x,names(x),c("seqnames","start","end"))
  }
  x[,GRanges(seqnames = seqnames,
             ranges = IRanges(
               start = start,
               end = end ))]
}

.rle_summit <- function(x)which.max(as.vector(x))

##' readBedFile
##'
##' Reads an increased bed file and adds the output into a \code{GRanges}
##' object.
##'
##' @param file a character with the increased bed file location.
##'
##' @export
##'
##' @return A \code{GRanges} object.
##'
##' @rdname readBedFile
##' @name readBedFile
##'
##' @references \url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}
##'
##' @examples
##'
##' reg = list.files(dr,pattern = "narrow",full.names =TRUE)
##' readBedFile(reg[1])
##'
##'
readBedFile <- function(file)
{
  dt = fread(file)
  gr = .dt2gr(dt)
  return(gr)
}

.data.table.GRanges <- function(x)
{
  dt <- data.table(seqnames = as.character( seqnames(x)),
    start = start(x),end = end(x),
    strand = as.character(strand(x)))
  return(dt)
}

.GRanges.data.table <- function(x)
{
  stopifnot(c("seqnames","start","end","strand") %in% names(x))
  return(GRanges(seqnames = x[,(seqnames)],
    ranges = .IRanges.data.table(x),
    strand = x[,(strand)]))
}

.IRanges.data.table <- function(x)
{
  stopifnot(c("seqnames","start","end","strand") %in% names(x))
  return(IRanges(start = x[,(start)],end = x[,(end)]))
}
#
# localMinima <- function(x)
# {
#   # Use -Inf instead if x is numeric (non-integer)
#   y <- diff(c(.Machine$integer.max, x)) > 0L
#   y <- cumsum(rle(y)$lengths)
#   y <- y[seq.int(1L, length(y), 2L)]
#   if (x[[1]] == x[[2]]) {
#     y <- y[-1]
#   }
#   return(y)
# }
#
# localMaxima <- function(x) localMinima(-x)
