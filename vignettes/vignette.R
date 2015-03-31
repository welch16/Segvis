## ----<style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------------------
  BiocStyle::latex()

## ----load,eval=FALSE--------------------------------------------------------------------
#  
#    library(Segvis)
#  

## ----real_load,include=FALSE,echo=FALSE,eval=TRUE---------------------------------------

  library(devtools)
  load_all("../")


## ----peaks_file,include=TRUE,echo=TRUE,eval=TRUE----------------------------------------

  peaks_file = "../inst/extdata/peaks/encode_K562_Ctcf_peaks_first3chr.narrowPeak"
  ctcf_peaks = read.table(peaks_file)
  head(ctcf_peaks,15)


## ----gr_peaks,include=TRUE,echo=TRUE,eval=TRUE------------------------------------------

  K = 2000
  ctcf_gr = GRanges(seqnames = ctcf_peaks$V1,
    ranges = IRanges(start = ctcf_peaks$V2,
      end = ctcf_peaks$V3),strand = "*")
  ctcf_gr = ctcf_gr[order(ctcf_peaks$V7,decreasing=TRUE)[1:K]]
  ctcf_gr


## ----parameters_segvis , include=TRUE,echo=TRUE,eval=TRUE, warning=FALSE----------------
  
  ctcf = Segvis(name = "ctcf_peaks",
    file = "../inst/extdata/reads/encode_K562_Ctcf_first3chr_Rep1.sort.bam",
    maxBandwidth = 101,fragLen = 200,isPET = FALSE,
    chr = c("chr1","chr2","chr3"))                 
  regions(ctcf) = ctcf_gr
  ctcf                 


## ----segvis_block , include=TRUE,echo=TRUE,eval=TRUE,message=FALSE----------------------
 
  ctcf = loadReads(ctcf, mc = 24)
  ctcf = matchReads(ctcf,mc = 24)
  ctcf = getCoverage(ctcf,mc = 24)
  ctcf_block = Segvis_block(ctcf,bw = 1,mc = 24)
  

## ----countReads,include=TRUE,echo=TRUE,eval=TRUE----------------------------------------

  ctcf_reads = countReads(ctcf)
  ctcf_reads
  

## ----normalize,include=TRUE,echo=TRUE,eval=TRUE-----------------------------------------

  normConst(ctcf_block) = ctcf_reads
  ctcf_block = normalize(ctcf_block)
  

## ----sessionInfo,include=TRUE,echo =TRUE,eval=TRUE,results="asis"-----------------------
  toLatex(sessionInfo())

