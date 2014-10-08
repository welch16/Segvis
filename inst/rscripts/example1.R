
rm(list = ls())
library(profile)

filedir = "../inst/extdata"
datadir = "../data"
rdatadir = "../inst/auxRData"

# How to build a GRanges object from narrowPeak format
bed_content = read.table(file = file.path(filedir,"encode_K562_Ctcf_peaks_first3chr.narrowPeak"),stringsAsFactors=FALSE)

peaks_ctcf = GRanges(seqnames = bed_content[,1],ranges =IRanges(start =bed_content[,2],end =bed_content[,3]),strand = '*')
peaks_ctcf$summit = start(peaks_ctcf) + bed_content[,10]

save(list = "peaks_ctcf",file = file.path(datadir,"peaks_ctcf.RData"))

# Initial parameters
name = "H3k27ac"
file = file.path(filedir,"encode_K562_H3k27ac_first3chr.sort.bam")
maxBw = 501
fl = 200

# Create profile
ourProfile = Profile(regionName = name,file = file,fileFormat = "bam",maxBandwidth = maxBw,fragLen = fl,remChr = "")

load(file = file.path(datadir,"peaks_ctcf.RData"))
windowExt = 1000
start(peaks_ctcf) = peaks_ctcf$summit - windowExt
end(peaks_ctcf) = peaks_ctcf$summit + windowExt
regions(ourProfile) = peaks_ctcf

mc =8
ourProfile = loadReads(ourProfile,mc)
ourProfile = matchReads(ourProfile,mc)
ourProfile = getCoverage(ourProfile,mc)

ourProfileMatrix1 = ProfileMatrix(ourProfile,1,mc)
ourProfileMatrix2 = ProfileMatrix(ourProfile,51,mc)
ourProfileMatrix3 = ProfileMatrix(ourProfile,301,mc)
ourProfileList = ProfileMatrixList(ourProfileMatrix1,ourProfileMatrix2,ourProfileMatrix3)
names(ourProfileList) = c("1","51","301")

p1 = plot.profiles(ourProfileList)
p2 = plot.profiles(ourProfileList,coord= -windowExt:windowExt)


save(list = c("p1","p2"),file = file.path(rdatadir,"bw_figs.RData"))






