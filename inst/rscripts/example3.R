
rm(list = ls())
library(profile)

filedir = "../inst/extdata"
datadir = "../data"
rdatadir = "../inst/auxRData"


load(file = file.path(datadir,"peaks_ctcf.RData"))
windowExt = 1000
start(peaks_ctcf) = peaks_ctcf$summit - windowExt
end(peaks_ctcf) = peaks_ctcf$summit + windowExt



names = c("H3k27ac","H3k4me1","H3k4me3")
files = list.files(filedir)
files_histones= file.path(filedir,do.call(c,lapply(names,function(x)files[grep(x,files)][2])))
maxBw = 501
fl = 200
mc =8

files = files[grep("Ctcf",files)]
files = files[grep(".sort.bam",files)]
files = files[!grepl(".bai",files)]
files = c(file.path(filedir,files[1]),files_histones)
names = c("Ctcf",names)

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

ourMatrices = lapply(ourMatrices,normalize.matrix)
ourList = ProfileMatrixList(ourMatrices)

bed_content = read.table(file = file.path(filedir,"encode_K562_dnase_Uw_1_first3chr.narrowPeak"),stringsAsFactors=FALSE)
dnase_regions = GRanges(seqnames = bed_content[,1],ranges =IRanges(start =bed_content[,2],end =bed_content[,3]),strand = '*')

bed_content = read.table(file = file.path(filedir,"encode_K562_Pol2b_first3chr.narrowPeak"),stringsAsFactors=FALSE)
pol2b = GRanges(seqnames = bed_content[,1],ranges =IRanges(start =bed_content[,2],end =bed_content[,3]),strand = '*')

ourMatrices = lapply(ourMatrices,function(x,dnase_regions,pol2b){
  regions(x)$dnase = countOverlaps(regions(x),dnase_regions)
  regions(x)$pol2b = countOverlaps(regions(x),pol2b)
  return(x)},dnase_regions,pol2b)

## save(list = "ourMatrices",file = "mm.RData")
## load("mm.RData")

p = plot.profiles(ourMatrices,coord = -windowExt:windowExt)
p1 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase > 0 & pol2b > 0)
p2 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase > 0 & pol2b == 0)
p3 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase == 0 & pol2b > 0)
p4 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase == 0 & pol2b == 0)

save(list  = c("p","p1","p2","p3","p4"),file = file.path(rdatadir,"subset_figs.RData"))
