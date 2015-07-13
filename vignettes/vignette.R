## ----<style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------------------
  BiocStyle::latex()

## ----load,eval=FALSE--------------------------------------------------------------------
#  
#    library(Segvis)
#  

## ----real_load,include=FALSE,echo=FALSE,eval=TRUE---------------------------------------

  ## this is to avoid all messages that appear when loading
  devtools::load_all("../")
  library(ggplot2)
  library(rbamtools)    
  library(Rsamtools)



## ----peaks_file,include=TRUE,echo=TRUE,eval=TRUE----------------------------------------

  peaks_file <- "../inst/extdata/peaks/encode_K562_Ctcf_peaks_first3chr.narrowPeak"
  ctcf_peaks <- read.table(peaks_file)
  head(ctcf_peaks,15)


## ----gr_peaks,include=TRUE,echo=TRUE,eval=TRUE------------------------------------------

  K <- 2000
  ctcf_gr <- GRanges(seqnames = ctcf_peaks$V1,
    ranges = IRanges(start = ctcf_peaks$V2,
      end = ctcf_peaks$V3),strand = "*")
  ctcf_gr <- ctcf_gr[order(ctcf_peaks$V7,decreasing=TRUE)[1:K]]
  ctcf_gr


## ----parameters_segvis , include=TRUE,echo=TRUE,eval=TRUE, warning=FALSE----------------
  
  ctcf <- buildSegvis(name = "ctcf_peaks",
    file = "../inst/extdata/reads/encode_K562_Ctcf_first3chr_Rep1.sort.bam",
    maxBandwidth = 101,fragLen = 200,isPET = FALSE,
    chr = c("chr1","chr2","chr3"))                 
  regions(ctcf) <- ctcf_gr
  ctcf                 


## ----segvis_block , include=TRUE,echo=TRUE,eval=TRUE,message=FALSE----------------------
 
  ctcf <- loadReads(ctcf, mc = 24)
  ctcf <- matchReads(ctcf,mc = 24)
  ctcf <- getCoverage(ctcf,mc = 24)
  ctcf_block <- Segvis_block(ctcf,bw = 1,mc = 24)
  

## ----countReads,include=TRUE,echo=TRUE,eval=TRUE----------------------------------------

  ctcf_reads <- countReads(ctcf)
  ctcf_reads
  

## ----normalize,include=TRUE,echo=TRUE,eval=TRUE-----------------------------------------

  normConst(ctcf_block) <- ctcf_reads
  ctcf_block <- normalize(ctcf_block)
  

## ----new_marks,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE-------------

  h3k27ac <- buildSegvis(name = "h3k27ac",
    file = "../inst/extdata/reads/encode_K562_H3k27ac_first3chr.sort.bam",
    maxBandwidth = 101,fragLen = 200,isPET = FALSE,
    chr = c("chr1","chr2","chr3"))
  regions(h3k27ac) <- ctcf_gr

  h3k27ac <- loadReads(h3k27ac, mc = 24)
  h3k27ac <- matchReads(h3k27ac,mc = 24)
  h3k27ac <- getCoverage(h3k27ac,mc = 24)
  h3k27ac_block <- Segvis_block(h3k27ac,bw = 1,mc = 24)

  h3k4me1 <- buildSegvis(name = "h3k4me1",
    file = "../inst/extdata/reads/encode_K562_H3k4me1_first3chr.sort.bam",
    maxBandwidth = 101,fragLen = 200,isPET = FALSE,
    chr = c("chr1","chr2","chr3"))
  regions(h3k4me1) <- ctcf_gr

  h3k4me1 <- loadReads(h3k4me1, mc = 24)
  h3k4me1 <- matchReads(h3k4me1,mc = 24)
  h3k4me1 <- getCoverage(h3k4me1,mc = 24)
  h3k4me1_block <- Segvis_block(h3k4me1,bw = 1,mc = 24)
                 

## ----block_list,include=TRUE,echo=TRUE,eval=TRUE,message =FALSE,warning=FALSE-----------

  h3k27ac_reads <- countReads(h3k27ac)  
  normConst(h3k27ac_block) <- h3k27ac_reads
  h3k27ac_block <- normalize(h3k27ac_block)

  h3k4me1_reads <- countReads(h3k4me1)                 
  normConst(h3k4me1_block) <- h3k4me1_reads
  h3k4me1_block <- normalize(h3k4me1_block)

  block_list <- Segvis_block_list(ctcf_block,h3k27ac_block,h3k4me1_block)
  names(block_list) <- c("ctcf","h3k27ac","h3k4me1")
                                   

## ----ex1_code,include=TRUE,echo=TRUE,eval=TRUE------------------------------------------

  rstart <- start(ctcf_gr)[1]
  rend <- end(ctcf_gr)[1]
  chr <- as.character(seqnames(ctcf_gr)[1])
  iden <- function(x)x  
  
  p1 <- plot_profiles(block_list,condition = seqnames == chr & start == rstart,
    coord = rstart:rend,FUN = iden,mc=24)
                 

## ----ex1_code2,include=TRUE,eval=TRUE,echo=TRUE,warning=FALSE,message=FALSE-------------

  p2 <- p1 + facet_grid(condition~.,scales = "free_y")
  p3 <- p2 + scale_colour_brewer(palette = "Dark2")+theme(legend.position = "none")


## ----ex1,include=TRUE,eval=TRUE,echo=FALSE,out.width='4.6cm',out.height='4cm',fig.show='hold'----

  p1
  p2
  p3            


## ----summits,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE-----------------------------

  summits <- findSummit(ctcf,bw=1,mc=24)
  ctcf_block <- addColumn(ctcf_block,name="summit",col=summits)
  ctcf_block
                    

## ----ex2_code,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE--------------

  window_ext <- 500
  new_start <- summits - window_ext
  new_end <- summits + window_ext
  new_regions <- GRanges(seqnames = seqnames(ctcf_gr),
    ranges = IRanges(start = new_start,end = new_end),strand ="*")
  str(width(new_regions))            

  regions(ctcf) <- new_regions
  regions(h3k27ac) <- new_regions
  regions(h3k4me1) <- new_regions
  all_segvis <- list("ctcf"=ctcf,"h3k27ac"=h3k27ac,"h3k4me1"=h3k4me1)
            

## ----new_profiles,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE----------

  do_all <- function(segvis_obj,bw,mc)
  {
    segvis_obj <- loadReads(segvis_obj,mc = mc)
    segvis_obj <- matchReads(segvis_obj,mc = mc)
    segvis_obj <- getCoverage(segvis_obj,mc = mc)
    out <- Segvis_block(segvis_obj,bw = bw , mc = mc)
    return(out)     
  }

  all_segvis_blocks <- lapply(all_segvis,do_all,bw = 1, mc = 24)
  nreads <- c(ctcf_reads,h3k27ac_reads,h3k4me1_reads)

  assign_and_normalize <- function(segvis_bl_obj,nreads)
  {
    normConst(segvis_bl_obj) <- nreads
    segvis_bl_obj <- normalize(segvis_bl_obj)
    return(segvis_bl_obj)    
  }
            
  all_segvis_blocks <- mapply(assign_and_normalize,
    all_segvis_blocks,nreads,SIMPLIFY=FALSE)
  all_segvis_blocks <- Segvis_block_list(all_segvis_blocks)
  names(all_segvis_blocks) <- names(all_segvis)
  

## ----ex2_plots,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE-------------

  q1 <- plot_profiles(all_segvis_blocks,FUN = mean,mc = 24,
    coord = -window_ext:window_ext)+xlab("distance to summit")+
    ylab("mean normalized counts")+
    scale_color_brewer(guide = guide_legend(title = "condition"),palette = "Dark2")+
    theme(legend.position = "top")+geom_vline(xintercept=0,linetype= 2)
            
  q2 <- plot_profiles(all_segvis_blocks,FUN = median,mc = 24,
    coord = -window_ext:window_ext)+xlab("distance to summit")+
    ylab("median normalized counts")+
    scale_color_brewer(guide = guide_legend(title = "condition"),palette = "Dark2")+
    theme(legend.position = "top")+geom_vline(xintercept=0,linetype= 2)    
                  

## ----ex2_plots2,include=TRUE,echo = TRUE,eval = TRUE,warning = FALSE,message = FALSE----

  varlog <- function(x)var(log(1 + x))
            
  q3 <- plot_profiles(all_segvis_blocks,FUN = varlog,mc = 24,
    coord = -window_ext:window_ext)+xlab("distance to summit")+
    ylab("variance of log( 1 + normalized counts)")+
    scale_color_brewer(guide = guide_legend(title = "condition"),palette = "Dark2")+
    theme(legend.position = "top")+geom_vline(xintercept=0,linetype= 2)


## ----ex2,include=TRUE,eval=TRUE,echo=FALSE,out.width='4.6cm',out.height='4cm',fig.show='hold'----

  q1
  q2
  q3            


## ----ex3_code,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE--------------

  # x is the genomic coordinates or distance to summit in thise case
  # y is the normalized counts
  # 1001 x 3 conditions = 3003 rows (in table below)
  new_data1 <- plot_data(all_segvis_blocks,FUN = mean,trim = .1,mc = 24,
    coord = -window_ext:window_ext)
  new_data2 <- plot_data(all_segvis_blocks,FUN = median,mc = 24,
    coord = -window_ext:window_ext)            
  new_data1
  

## ----ex3_plots,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,message=FALSE-------------

  p4 <- p3 %+% new_data1 + ylab("average normalized counts")
  p5 <- p3 %+% new_data2 + ylab("median normalized counts")
             

## ----ex3,include=TRUE,eval=TRUE,echo=FALSE,out.width='4.6cm',out.height='4cm',fig.show='hold'----

  p4
  p5            


## ----ex4_code,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE--------------

  dnase_file <- "../inst/extdata/peaks/encode_K562_dnase_openChrom_first3chr.narrowPeak"
  list.files("../inst/extdata/peaks/")
  dnase_sites = read.table(dnase_file)
  dnase_gr <- GRanges(seqname = dnase_sites$V1,
    ranges = IRanges(start = dnase_sites$V2,end = dnase_sites$V3),
    strand = "*")
 dnase_gr


## ----ex4_code2,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE-------------

  nr_overlaps <- countOverlaps(regions(all_segvis_blocks[[1]]),dnase_gr)
  all_segvis_blocks <- lapply(all_segvis_blocks,
    addColumn,name = "dnase_overlaps",col =nr_overlaps)
  all_segvis_blocks[[1]]
                          

## ----ex4_code3,subset_example,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE----

  ctcf_subset <- subset(all_segvis_blocks[[1]],dnase_overlaps > 0)
  ctcf_subset
  cover_table(ctcf_subset)
  

## ----ex4_plots,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE-------------

  
  s1 <- plot_profiles(all_segvis_blocks,FUN = mean,mc = 24,
    condition = dnase_overlaps > 0 ,    
    coord = -window_ext:window_ext)+xlab("distance to summit")+
    ylab("mean normalized counts")+
    scale_color_brewer(guide = guide_legend(title = "condition"),palette = "Dark2")+
    theme(legend.position = "top")+geom_vline(xintercept=0,linetype= 2)
            
  s2 <- plot_profiles(all_segvis_blocks,FUN = median,mc = 24,
    condition = dnase_overlaps > 0 ,    
    coord = -window_ext:window_ext)+xlab("distance to summit")+
    ylab("median normalized counts")+
    scale_color_brewer(guide = guide_legend(title = "condition"),palette = "Dark2")+
    theme(legend.position = "top")+geom_vline(xintercept=0,linetype= 2)

  s3 <- plot_profiles(all_segvis_blocks,FUN = varlog,mc = 24,
    condition = dnase_overlaps > 0 ,
    coord = -window_ext:window_ext)+xlab("distance to summit")+
    ylab("variance of log( 1 + normalized counts)")+
    scale_color_brewer(guide = guide_legend(title = "condition"),palette = "Dark2")+
    theme(legend.position = "top")+geom_vline(xintercept=0,linetype= 2)
                       

## ----ex4,include=TRUE,eval=TRUE,echo=FALSE,out.width='4.6cm',out.height='4cm',fig.show='hold'----

  s1
  s2
  s3


## ----ex5_code,include=TRUE,echo=TRUE,eval=TRUE,message=FALSE,warning=FALSE--------------

  mean_overlap_data <- plot_data(all_segvis_blocks,FUN = mean,mc = 24,
    condition = dnase_overlaps > 0,
    coord = -window_ext:window_ext)
  mean_comp_data <- plot_data(all_segvis_blocks,FUN = mean,mc = 24,
    condition = dnase_overlaps ==  0,
    coord = -window_ext:window_ext)
  new_data <- rbind(mean_overlap_data[,overlap:="yes"],
    mean_comp_data[,overlap:="no"])
  fancy_plot <- ggplot(new_data,aes(x,y,colour = condition))+geom_line(size=1.1)+
    facet_grid(overlap~.,scales = "free_y")+theme(legend.position = "top")+
    ggtitle("DHS overlaps")+geom_vline(xintercept = 0,linetype=2,size=1.1)+
    xlab("distance to summit")+ylab("average coverage")+
    scale_color_brewer(palette = "Set1")
                            

## ----ex5,include=TRUE,eval=TRUE,echo=FALSE,out.width='6cm',out.height='5cm',fig.show='hold'----
  fancy_plot

## ----sessionInfo,include=TRUE,echo =TRUE,eval=TRUE,results="asis"-----------------------
  toLatex(sessionInfo())

