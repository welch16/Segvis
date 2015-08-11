 
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(Segvis)
library(parallel)

files_dir <- "inst/extdata"
mc <- detectCores()

peak_dir <- file.path(files_dir,"peaks")
read_dir <- file.path(files_dir,"reads")
example_dir <- file.path(files_dir,"example")


ctcf <- "encode_K562_Ctcf_peaks_first3chr.narrowPeak"
peaks <- read.table(file.path(peak_dir,ctcf))
peaks <- data.table(peaks)

peaks <- split(peaks,peaks[,(V1)])

set.seed(12345)

N <- 50
rows <- lapply(peaks,function(x,N)sample(nrow(x),N),N)

out_peaks <- do.call(rbind,
    mapply(function(peak,idx)peak[idx],peaks,rows,SIMPLIFY = FALSE))

write.table(out_peaks,file = file.path(example_dir,ctcf),quote = FALSE,
    row.names = FALSE,col.names = FALSE,sep = "\t")            

reads_files <- list.files(read_dir)
reads_files <- reads_files[grep("bai",reads_files,invert = TRUE)]
reads_files <- reads_files[grep("sort",reads_files)]

## create segvis object

ctcf <- buildSegvis(name = "ctcf_peaks",
  file = file.path(read_dir,reads_files[1]),
  maxBandwidth = 101,fragLen = 200,isPET = FALSE,
  chr = c("chr1","chr2","chr3"))
regions(ctcf) <- GRanges(seqnames = out_peaks[,(V1)],
    ranges = IRanges(start = out_peaks[,(V2)],
      end = out_peaks[,(V3)]),
    strand = "*")                         


ctcf <- loadReads(ctcf,mc = mc)
ctcf <- matchReads(ctcf,mc = mc)
ctcf <- getCoverage(ctcf,mc = mc)
ctcf_block <- Segvis_block(ctcf,bw = 1,mc = mc)
summits <- findSummit(ctcf,bw = 1,mc = mc)

window_ext <- 500
regions <- regions(ctcf)

centered <- regions
start(centered) <- summits - window_ext
end(centered) <- summits + window_ext

regions <- reduce(c(centered,regions))

OUT <- mapply(filterBam, file.path(read_dir,reads_files),file.path(example_dir,reads_files),
  MoreArgs = list(param = ScanBamParam(which = regions)),SIMPLIFY = FALSE)

peakfiles <- list.files(peak_dir)
peakfiles <- peakfiles[grep("Ctcf",peaks,invert = TRUE)]

peaks <- lapply(file.path(peak_dir,peakfiles),function(x)data.table(read.table(x)))

gr <- lapply(peaks,function(x){
    GRanges(seqnames = x[,(V1)],
        ranges = IRanges(start = x[,(V2)],end = x[,(V3)]),
        strand = "*")})

overlaps <- lapply(gr,findOverlaps,regions)

reduced_peaks <- mapply(function(peak,ov)peak[queryHits(ov)],peaks,overlaps,SIMPLIFY = FALSE)

u <- mapply(function(peak,file){
  write.table(peak,file = file,sep = "\t", row.names = FALSE, col.names = FALSE,
    quote = FALSE)},reduced_peaks,file.path(example_dir,peakfiles),SIMPLIFY = FALSE)



