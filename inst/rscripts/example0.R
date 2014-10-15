
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
maxBw = 1
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


ourMatrices = lapply(ourProfiles,ProfileMatrix,1,mc)
names(ourMatrices) = names



ourMatrices = lapply(ourMatrices,normalize.matrix)
ourList = ProfileMatrixList(ourMatrices)



save(file=file.path(rdatadir,"example0.RData"),list = "ourList")
  

starts = start(regions(ourList[[1]]))
ends = end(regions(ourList[[1]]))
seqs =as.character( seqnames(regions(ourList[[1]])))

i=10
p1 = plot.profiles(ourList,condition = seqnames == "chr1" & start == starts[i],coord = seq(starts[i],ends[i]))
p1



save(list  = c("p","p1","p2","p3","p4"),file = file.path(rdatadir,"subset_figs.RData"))
