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
  
  ctcf_gr = GRanges(seqnames = ctcf_peaks$V1,
    ranges = IRanges(start = ctcf_peaks$V2,
      end = ctcf_peaks$V3),strand = "*")
  ctcf_gr


## ----parameters_segvis , include=TRUE,echo=TRUE,eval=TRUE, warning=FALSE----------------
  
  ctcf = Segvis(name = "ctcf_peaks",
    file = "../inst/extdata/reads/encode_K562_Ctcf_first3chr_Rep1.bam",
    maxBandwidth = 101,fragLen = 200,isPET = FALSE,
    chr = c("chr1","chr2","chr3"))                 
  regions(ctcf) = ctcf_gr
  ctcf                 


## ----sessionInfo,include=TRUE,echo =TRUE,eval=TRUE,results="asis"-----------------------
  toLatex(sessionInfo())

