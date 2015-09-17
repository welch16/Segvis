##' @import data.table
##' @import rbamtools
##' @import Rsamtools
##' @importFrom GenomicAlignments readGAlignments
##' @importFrom GenomicAlignments seqnames
##' @importFrom GenomicAlignments findOverlaps
##' @importFrom GenomicAlignments start
##' @importFrom GenomicAlignments end
##' @importFrom GenomicRanges GRanges
##' @importFrom GenomicRanges start<-
##' @importFrom GenomicRanges end<-
##' @importFrom GenomicRanges elementMetadata
##' @importFrom GenomicRanges elementMetadata<-
##' @importFrom GenomicRanges resize
##' @importFrom GenomicRanges coverage
##' @importFrom GenomicRanges queryHits
##' @importFrom GenomicRanges subjectHits
##' @importFrom GenomicRanges width
##' @importFrom GenomicRanges strand
##' @importFrom IRanges IRanges
##' @importFrom IRanges IRangesList
##' @importFrom GenomeInfoDb seqlengths
##' @importFrom S4Vectors nrun
##' @importFrom S4Vectors runValue
##' @importFrom S4Vectors runLength
##' @import parallel
NULL

##' @rdname name-methods
##' @aliases name
##' @docType methods
##' @exportMethod name
setMethod("name",
  signature = signature(object = "segvis"),
  definition = function(object)object@name
)           

##' @rdname regions-methods
##' @aliases regions
##' @docType methods
##' @exportMethod regions
setMethod("regions",
  signature = signature(object = "segvis"),
  definition = function(object)object@regions
)           

##' @rdname file-methods
##' @aliases file
##' @docType methods
##' @exportMethod file
setMethod("file",
  signature = signature(object = "segvis"),
  definition = function(object)object@file
)           

##' @rdname maxBandwidth-methods
##' @aliases maxBandwidth
##' @docType methods
##' @exportMethod maxBandwidth
setMethod("maxBandwidth",
  signature = signature(object = "segvis"),
  definition  = function(object)object@maxBandwidth
)           

##' @rdname fragLen-methods
##' @aliases fragLen
##' @docType methods
##' @exportMethod fragLen
setMethod("fragLen",
  signature = signature(object = "segvis"),
  definition = function(object)object@fragLen
)           

##' @rdname chr-methods
##' @aliases chr
##' docType methods
##' @exportMethod chr
setMethod("chr",
  signature = signature(object = "segvis"),
  definition = function(object)object@chr
)          

##' @rdname isPET-methods
##' @aliases isPET
##' @docType methods
##' @exportMethod isPET
setMethod("isPET",
  signature = signature(object = "segvis"),
  definition = function(object)object@isPET
)          

##' @rdname readsF-methods
##' @docType methods
##' @aliases readsF
##' @exportMethod readsF
setMethod("readsF",
  signature = signature(object = "segvis"),
  definition = function(object)readsF(object@reads)
)           


##' @rdname readsR-methods
##' @docType methods
##' @aliases readsR
##' @exportMethod readsR
setMethod("readsR",
  signature = signature(object = "segvis"),
  definition = function(object)readsR(object@reads)
)

##' @rdname profiles-methods
##' @docType methods
##' @aliases profiles
##' @exportMethod profiles
setMethod("profiles",
  signature = signature(object = "segvis"),
  definition = function(object)object@profiles
)  

## Set methods

##' @rdname name-methods
##' @aliases name<-
##' @docType methods
##' @exportMethod name<-
setReplaceMethod("name",
  signature = signature(object = "segvis",value = "character"),
  definition = function(object,value){
    object@name <- value
    return(object)
})
   

##' @rdname regions-methods
##' @aliases regions<-
##' @docType methods
##' @exportMethod regions<-
setReplaceMethod("regions",
  signature = signature(object = "segvis",value = "GRanges"),
  definition = function(object,value){
    
    stopifnot(class(value) == "GRanges")
    
    unique_chr <- unique(as.character(seqnames(value)))
    if(length(chr(object)) > length(unique_chr)){
      warning("There are more chromosomes in chr(object) than in the regions. User's supplied additional chromosomes are being removed")
      chr(object) <- chr(object)[chr(object) %in% unique_chr]
    }
    if(all(chr(object) != "")){
      chr_regions <- unique(as.character(seqnames(value)))
      chr_in <- sapply(chr_regions,function(x)x%in%chr(object))
      if(!all(chr_in)){
        warning("Removing reads that aren't in user's defined chr")
      }
      value_dt <- .data.table.GRanges(value)
      setkey(value_dt,seqnames)
      value <- .GRanges.data.table(value_dt[chr(object)])             
    }
    object@regions <- value
    object@.haveRegions <- TRUE
    return(object)
})

##' @rdname file-methods
##' @aliases file<-
##' @docType methods
##' @exportMethod file<-
setReplaceMethod("file",
  signature = signature(object = "segvis",value = "character"),
  definition = function(object,value){

    stopifnot(class(value) == "character")

    object@file <- value
    return(object)
})    



##' @rdname maxBandwidth-methods
##' @aliases maxBandwidth<-
##' @docType methods
##' @exportMethod maxBandwidth<-
setReplaceMethod("maxBandwidth",
  signature = signature(object = "segvis", value = "numeric"),
  definition = function(object,value){
    
    stopifnot(value >= 1)
    stopifnot(value %% 2 == 1)
    
    object@maxBandwidth <- value
    return(object)
})    

##' @rdname fragLen-methods
##' @aliases fragLen<-
##' @docType methods
##' @exportMethod fragLen<-
setReplaceMethod("fragLen",
  signature = signature(object = "segvis", value = "numeric"),
  definition = function(object,value){
    stopifnot(value >= 0)
    stopifnot(value == floor(value))
    object@fragLen <- value
    return(object)
})    

##' @rdname chr-methods
##' @aliases chr<-
##' @docType methods
##' @exportMethod chr<-
setReplaceMethod("chr",
  signature = signature(object = "segvis",value = "character"),
  definition = function(object,value){
    stopifnot(is.character(value))
    object@chr <- value
    return(object)
})

##' @rdname isPET-methods
##' @aliases isPET<-
##' @docType methods
##' @exportMethod isPET<-
setReplaceMethod("isPET",
  signature = signature(object = "segvis",value = "logical"),
  definition = function(object,value){
    stopifnot(is.logical(value))
    object@isPET <- value
    return(object)
})    

# @rdname methods-segvis-show
# @name show
setMethod("show",
  signature = signature(object = "segvis"),
  definition = function(object){
#    cat("---------------------------\n")
    cat("Segvis for",name(object),"regions\n")
    cat("Paired-end Tags: ",isPET(object),"\n")
    cat("Fragment length:",fragLen(object),"\n")
    cat("Max Bandwidth:", maxBandwidth(object),"\n")
    cat("Using reads file:\n")
    cat(file(object),sep = "\n")
    len <- length(seqlengths(regions(object)))
    if( len > 0){
      cat("Using regions for",len,"chromosomes\n")
      show(regions(object))      
    }else{
      cat("**Not regions loaded**\n")
    }
#    cat("---------------------------\n")
})


##' @rdname readsF-methods
##' @docType methods
##' @aliases readsF<-
##' @exportMethod readsF
setReplaceMethod("readsF",
  signature = signature(object = "segvis",value = "list"),
  definition = function(object,value){
    readsF(object@reads) = value
    return(object)
  }
)                 


##' @rdname readsR-methods
##' @docType methods
##' @aliases readsR<-
##' @exportMethod readsR
setReplaceMethod("readsR",
  signature = signature(object = "segvis",value = "list"),
  definition = function(object,value){
    readsR(object@reads) <- value
    return(object)
  }
)                 

##' @rdname loadReads-methods
##' @aliases loadReads
##' @docType methods
##' @exportMethod loadReads
setMethod("loadReads",
  signature = signature(object = "segvis",mc = "numeric"),
  definition = function(object,mc ){
    ## reads the bam file
    message("Reading ",file(object))

    bai <- paste0(file(object),".bai")
    if(!file.exists(bai)){
      warning("Creating index file ",bai)
      message("Creating index file ",bai)
      indexBam(file(object))
    }
    
    side <- (maxBandwidth(object)-1)/2      
    regions_to_load <- regions(object)    
    start(regions_to_load) <- start(regions_to_load) - side - fragLen(object)
    end(regions_to_load) <- end(regions_to_load) + side + fragLen(object)          
    
    if(isPET(object)){
      message("Setting PET flag")
      ## when the reads are paired end tags, it gets the qname to match the
      ## both ends of the fragment
      pet_flag <- scanBamFlag(isPaired = TRUE)
      param <- ScanBamParam(which = regions_to_load,flag = pet_flag,what = "qname")
    }else{
      param <- ScanBamParam(which = regions_to_load)
    }    
    greads <- readGAlignments(file(object), param = param,use.names = FALSE)
    if(length(greads) == 0){
      warning("Can't read ",file(object), " considering regions, going to try without them")
      param@which <- IRangesList()
      greads <- readGAlignments(file(object), param = param,use.names = FALSE)
      overlap <- findOverlaps(as(greads,"GRanges"),regions_to_load)
      greads <- greads[queryHits(overlap)]
      rm(overlap)      
    }
    if(isPET(object)){
      ## convert the qname into a numeric value for computation efficiency
      qname <- as.numeric(as.factor(elementMetadata(greads)[["qname"]]))
      greads <- .data.table.GRanges(as(greads, "GRanges"))
      greads[,name:=qname] # add qname to greads
      setorder(greads,name)
    }else{
      greads <- .data.table.GRanges(as(greads, "GRanges"))
    }
    setkey(greads,seqnames)    
    message("Bam file loaded")

    ## validates the chromosomes, that must coincide between
    ## regions, reads and the ones supplied by user
    chr_reads <- unique(greads$seqnames)
    chr_in <- sapply(chr_reads,function(x)x%in%chr(object))
    if(!all(chr_in)){
      warning("Removing reads that aren't in user's supplied chromosome")
    }

    ## update options of data.table for when there are no fragment for a region,
    ## then it would get an empty row for that chromosome (other value is NA)

    chr_in <- sapply(chr(object),function(x)x%in%chr_reads)
    if(!all(chr_in)){
      missing_chr <- names(which(!chr_in))
      for(ch in missing_chr){
        warning("There are no reads for ",ch," in " ,file(object))
      }        
    }
    
    options(datatable.nomatch =0)
    greads <- greads[chr(object)]

    ## separates by chromosome and strand
    message("Separating reads by strand")
    message("Forward strand reads extracted")
    greads1 <- separate.by.chrom(greads,chr(object),"+",mc,sort=TRUE)
    names(greads1) <- chr(object)
    message("Reverse strand reads extracted")
    greads2 <- separate.by.chrom(greads,chr(object),"-",mc,sort=TRUE)
    names(greads2) <- chr(object)
    message("Finished separating reads")

    ## Creates the reads object
    object@reads <- new("reads",readsF = greads1,readsR = greads2)
    object@.haveReads <- TRUE
    message("Reading bam files... Done")
    return(object)
})

##' @rdname matchReads-methods
##' @aliases matchReads
##' @docType methods
##' @exportMethod matchReads
setMethod("matchReads",
  signature = signature(object = "segvis",mc = "numeric"),
  definition = function(object,mc = 8){
    if(object@.haveReads & object@.haveRegions){

      ## separates the regions by chromosome, extend the regions by (maxBw - 1)/2 to each
      ## side
      
      side <- (maxBandwidth(object)-1)/2      
      chr <- names(seqlengths(regions(object)))
      match_regions <- regions(object)
      start(match_regions) <- start(match_regions) - side - fragLen(object)
      end(match_regions) <- end(match_regions) + side + fragLen(object)
      match_regions <- separate.by.chrom(.data.table.GRanges(match_regions),
        chr, "*",mc,sort=FALSE)
      names(match_regions) <- chr

      ## find overlaps between reads and regions   
      overlapsF <- mcmapply(.find.overlaps,readsF(object),
        match_regions,SIMPLIFY=FALSE,mc.cores =mc,mc.silent =TRUE,mc.preschedule=TRUE)
      overlapsR <- mcmapply(.find.overlaps,readsR(object),
        match_regions,SIMPLIFY=FALSE,mc.cores =mc,mc.silent =TRUE,mc.preschedule=TRUE)
              
      ## matching reads
      message("Matching reads for forward strand")
      readsF(object) <- mcmapply(.match.reads,readsF(object),overlapsF,
        SIMPLIFY= FALSE,mc.cores = mc ,mc.silent =TRUE,mc.preschedule=TRUE)
      message("Matching reads for reverse strand")
      readsR(object) <- mcmapply(.match.reads,readsR(object),overlapsR,
        SIMPLIFY= FALSE,mc.cores = mc ,mc.silent =TRUE,mc.preschedule=TRUE)

      object@.readsMatched <- TRUE
     return(object)      
    }else{
      warning("Check that both reads and regions are loaded")
    }
})

##' @rdname getCoverage-methods
##' @aliases getCoverage
##' @docType methods
##' @exportMethod getCoverage
setMethod("getCoverage",
  signature = signature(object = "segvis",mc = "numeric"),
  definition = function(object, mc = 8){
    if(object@.readsMatched == TRUE){
      # init coverage calculation
      chr <- names(seqlengths(regions(object)))
      match_regions <- separate.by.chrom(.data.table.GRanges(regions(object)),
          chr, "*",mc,sort=FALSE)
      nregions <- mclapply(match_regions,nrow,mc.cores = mc,mc.silent =TRUE)

      # coverage calculation
      message("Calculating coverage")      
      curves <- mapply(calculate_chrom_coverage,chr,nregions,
        MoreArgs = list(object,mc),SIMPLIFY=FALSE)      
      names(curves) <- chr
      object@profiles <- curves
      object@.coverageCalculated <- TRUE
      message("Coverage done")
      return(object)                
    }else{
      warning("The reads haven't been matched yet")
    }
})    

##' @rdname findSummit-methods
##' @aliases findSummit
##' @docType methods
##' @exportMethod findSummit
setMethod("findSummit",          
  signature = signature(object = "segvis",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){

    message("Calculating profiles for each chromosome")
    profile_curves <- .calculate_profile_curves(object,bw,mc)

    message("Finding summits for each chromosome")
    chr <- names(seqlengths(regions(object)))
    match_regions <- separate.by.chrom(.data.table.GRanges(regions(object)),
        chr, "*",mc,sort=FALSE)    
    names(match_regions) <- chr
    
    ## find summits for each chromosome
    summits_chr <- lapply(chr,function(chrom,regions,curves){
      message("Finding summit for regions of ",chrom)
      reg <- regions[[chrom]]
      regionStart <- reg[,(start)]
      regionEnd <- reg[,(end)]
      curve <- curves[[chrom]]
      summits <- mcmapply(.find_summit,curve,regionStart,regionEnd,
        SIMPLIFY=FALSE,mc.cores = mc,mc.silent = TRUE,mc.preschedule=TRUE)
      return(unlist(summits))         
  },match_regions,profile_curves)
  return(unlist(summits_chr))    
})

##' @rdname countReads-methods
##' @aliases countReads
##' @docType methods
##' @exportMethod countReads
setMethod("countReads",
  signature = signature(object = "segvis"),
  definition = function(object){
    ## Check the case for a PET file
    reader <- bamReader(file(object),idx=TRUE)
    readMatrix <- bamCountAll(reader)
    counts <- sum(readMatrix$nAligns)
    if(isPET(object)){
      counts <- ceiling(counts/2)
    }
    return(counts)  
})

##' @rdname joinProfiles-methods
##' @aliases joinProfiles
##' @docType methods
##' @exportMethod joinProfiles
setMethod("joinProfiles",
  signature = signature(object = "segvis",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){

  ## calculate the coverage vectors
  message("Calculating profiles for each chromosome")
  profile_curves <- .calculate_profile_curves(object,bw,mc)

  ## init regions data
  chr <- names(seqlengths(regions(object)))
  match_regions <- separate.by.chrom(.data.table.GRanges(regions(object)),
    chr, "*",mc,sort=FALSE)    
  names(match_regions) <- chr

  ## join regions a profiles into data.table's
  joined_info <- mapply(.join_info,chr,match_regions,profile_curves,
    MoreArgs <- list(mc),SIMPLIFY=FALSE)
  joined_info <- do.call(rbind,joined_info)
  
  message("All profiles joint")
  return(joined_info)
})    

##' @rdname Segvis_block-methods
##' @aliases Segvis_block
##' @docType methods
##' @exportMethod Segvis_block
setMethod("Segvis_block",
  signature = signature(object = "segvis",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc){
    cover_table <- joinProfiles(object,bw,mc)    
    nm <- name(object)
    gr <- regions(object)
    segvis_block <- new("segvis_block",name = nm,regions = gr,
      cover_table = cover_table,bandwidth = bw,normConst = 1)
    return(segvis_block)
})             
