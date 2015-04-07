
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

region_coverage <-  function(fwd,bwd,isPET,fragLen,chr)
{
  if(isPET){
    qnames = intersect(fwd[,(name)],bwd[,(name)])
    setkey(fwd,name)
    setkey(bwd,name)
    pair_fwd = fwd[name %in% qnames]
    pair_bwd = bwd[name %in% qnames]
    if(!identical(pair_fwd[,(name)],pair_bwd[,(name)])){
      setorder(pair_fwd,name)
      setorder(pair_bwd,name)
    }
    reads = copy(pair_fwd)
    reads[,end:= pair_bwd[,(end)]]
    reads[,name:=NULL]
    reads[,strand:="*"]
    reads = reads[end - start > -1]  
    reads = .GRanges.data.table(reads)
  }else{
    if(fragLen == 0){
      reads = .GRanges.data.table(rbind(fwd,bwd))
    }else{
      reads = resize(.GRanges.data.table(rbind(fwd,bwd)),fragLen)
    }
  }
  return(coverage(reads)[[chr]])
}


calculate_chrom_coverage <- function(chr,nreg,object,mc)
{
  fwd_reads = readsF(object)[[chr]]
  bwd_reads = readsR(object)[[chr]]
  setkey(fwd_reads,match)
  setkey(bwd_reads,match)
  sep_fwd_reads = mclapply(1:nreg,function(i,reads)reads[match==i,],
    fwd_reads,mc.cores = mc,mc.silent= TRUE)
  sep_bwd_reads = mclapply(1:nreg,function(i,reads)reads[match==i,],
    bwd_reads,mc.cores = mc,mc.silent= TRUE)
  curves = mcmapply(region_coverage,sep_fwd_reads,sep_bwd_reads,
    MoreArgs = list(isPET(object),fragLen(object),chr),SIMPLIFY =FALSE,
    mc.cores = mc,mc.silent = TRUE)

  return(curves)
}

.join_info <- function(chr,region,profile,mc)
{
  nreg = nrow(region)
  chr_cover = mclapply(1:nreg,function(i,region,profile)
    data.table(chr = chr,match=i,coord = seq(region[i,(start)],region[i,(end)],by=1),
      tagCounts = profile[[i]]),region,profile,mc.cores = mc,mc.silent = mc)
  chr_cover = do.call(rbind,chr_cover)
  return(chr_cover)
}
