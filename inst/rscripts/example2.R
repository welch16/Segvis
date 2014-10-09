
rm(list = ls())
library(profile)

filedir = "../inst/extdata"
datadir = "../data"
rdatadir = "../inst/auxRData"

# How to build a GRanges object from narrowPeak format
bed_content = read.table(file = file.path(filedir,"encode_K562_Ctcf_peaks_first3chr.narrowPeak"),stringsAsFactors=FALSE)
peaks_ctcf = GRanges(seqnames = bed_content[,1],ranges =IRanges(start =bed_content[,2],end =bed_content[,3]),strand = '*')
peaks_ctcf$summit = start(peaks_ctcf) + bed_content[,10]

names = c("H3k27ac","H3k4me1","H3k4me3")
files = list.files(filedir)
files = file.path(filedir,do.call(c,lapply(names,function(x)files[grep(x,files)][2])))
maxBw = 501
fl = 200
mc =8

load(file = file.path(datadir,"peaks_ctcf.RData"))
windowExt = 1000
start(peaks_ctcf) = peaks_ctcf$summit - windowExt
end(peaks_ctcf) = peaks_ctcf$summit + windowExt
regions(ourProfile) = peaks_ctcf


  
ourProfiles = mapply(Profile,names,files,MoreArgs = list("bam",maxBw,fl,""),SIMPLIFY=FALSE )
names(ourProfiles) = names
ourProfiles = lapply(ourProfiles,function(x,peaks_ctcf){
  regions(x) = peaks_ctcf
  return(x)},peaks_ctcf)
ourProfiles = lapply(ourProfiles,loadReads,mc)
ourProfiles = lapply(ourProfiles,matchReads,mc)
ourProfiles = lapply(ourProfiles,getCoverage,mc)

ourMatrices = lapply(ourProfiles,ProfileMatrix,251,mc)
names(ourMatrices) = names
ourList_notScaled = ProfileMatrixList(ourMatrices)

q1 = plot.profiles(ourList_notScaled,coord= -windowExt:windowExt)

ourMatrices = lapply(ourMatrices,normalize.matrix)
ourList = ProfileMatrixList(ourMatrices)

q2 = plot.profiles(ourList,coord= -windowExt:windowExt)
q3 = plot.profiles(ourList,coord= -windowExt:windowExt,trim = 1)
q4 = plot.profiles(ourList,coord= -windowExt:windowExt,trim = .25)


save(list = c("q1","q2","q3","q4"),file = file.path(rdatadir,"trim_figs.RData"))



