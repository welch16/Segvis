
# Methods for profile class

## Get methods

#' @rdname methods-segvis-gs
#' @name name 
setMethod("name",
  signature = signature(object = "segvis"),
  definition = function(object)object@name
)           

#' @rdname methods-segvis-gs
#' @name regions
setMethod("regions",
  signature = signature(object = "segvis"),
  definition = function(object)object@regions
)           

#' @rdname methods-segvis-gs
#' @name file
setMethod("file",
  signature = signature(object = "segvis"),
  definition = function(object)object@file
)           

#' @rdname methods-segvis-gs
#' @name maxBandwidth
setMethod("maxBandwidth",
  signature = signature(object = "segvis"),
  definition  = function(object)object@maxBandwidth
)           

#' @rdname methods-segvis-gs
#' @name fragLen
setMethod("fragLen",
  signature = signature(object = "segvis"),
  definition = function(object)object@fragLen
)           

#' @rdname methods-segvis-gs
#' @name chr
setMethod("chr",
  signature = signature(object = "segvis"),
  definition = function(object)object@chr
)          

#' @rdname methods-segvis-gs
#' @name isPET
setMethod("isPET",
  signature = signature(object = "segvis"),
  definition = function(object)object@isPET
)          

#' @rdname methods-segvis-gs
#' @name readsF
setMethod("readsF",
  signature = signature(object = "segvis"),
  definition = function(object)readsF(object@reads)
)           

#' @rdname methods-segvis-gs
#' @name readsR
setMethod("readsR",
  signature = signature(object = "segvis"),
  definition = function(object)readsR(object@reads)
)

# @rdname methods-segvis-gs
# @name matchF
## setMethod("matchF",
##   signature = signature(object = "segvis"),
##   definition = function(object)matchF(object@match)
## )

# @rdname methods-segvis-gs
# @name matchR
## setMethod("matchR",
##   signature = signature(object = "segvis"),
##   definition = function(object)matchR(object@match)
## )

#' @rdname methods-segvis-gs
#' @name profileCurve
setMethod("profileCurve",
  signature = signature(object = "segvis"),
  definition = function(object)object@profileCurve
)  



## Set methods

#' @rdname methods-segvis-gs
#' @name name<-
setReplaceMethod("name",
  signature = signature(object = "segvis",value = "character"),
  definition = function(object,value){
    object@name = value
    return(object)
})
   
#' @rdname methods-segvis-gs
#' @name regions<-
setReplaceMethod("regions",
  signature = signature(object = "segvis",value = "GRanges"),
  definition = function(object,value){
    stopifnot(class(value) == "GRanges")
    unique_chr = unique(as.character(seqnames(value)))
    if(length(chr(object)) >
      length(unique_chr)){
      warning("There are more chromosomes in chr(object) than in the regions. User's supplied additional chromosomes are being removed")
      chr(object) = chr(object)[chr(object) %in% unique_chr]
    }
    if(all(chr(object) != "")){
      chr_regions = unique(as.character(seqnames(value)))
      chr_in =sapply(chr_regions,function(x)x%in%chr(object))
      if(!all(chr_in)){
        warning("Removing reads that aren't in user's defined chr")
      }
      value_dt = .data.table.GRanges(value)
      setkey(value_dt,seqnames)
      value = .GRanges.data.table(value_dt[chr(object)])             
    }
    object@regions = value
    object@.haveRegions = TRUE
    return(object)
})

#' @rdname methods-segvis-gs
#' @name file<-
setReplaceMethod("file",
  signature = signature(object = "segvis",value = "character"),
  definition = function(object,value){
    stopifnot(class(value) == "character")
    object@file = value
    return(object)
})    

#' @rdname methods-segvis-gs
#' @name maxBandwidth<-
setReplaceMethod("maxBandwidth",
  signature = signature(object = "segvis", value = "numeric"),
  definition = function(object,value){
    stopifnot(value >= 1)
    stopifnot(value %% 2 == 1)
    object@maxBandwidth = value
    return(object)
})    

#' @rdname methods-segvis-gs
#' @name fragLen<-
setReplaceMethod("fragLen",
  signature = signature(object = "segvis", value = "numeric"),
  definition = function(object,value){
    stopifnot(value >= 0)
    stopifnot(value == floor(value))
    object@fragLen = value
    return(object)
})    

#' @rdname methods-segvis-gs
#' @name chr<-
setReplaceMethod("chr",
  signature = signature(object = "segvis",value = "character"),
  definition = function(object,value){
    stopifnot(is.character(value))
    object@chr = value
    return(object)
})

#' @rdname methods-segvis-gs
#' @name isPET<-
setReplaceMethod("isPET",
  signature = signature(object = "segvis",value = "logical"),
  definition = function(object,value){
    stopifnot(is.logical(value))
    object@isPET = value
    return(object)
})    

# @rdname methods-segvis-show
# @name show
setMethods("show",
  signature = signature(object = "segvis"),
  definition = function(object){
#    cat("---------------------------\n")
    cat("Profile for",name(object),"regions\n")
    cat("Paired-end Tags: ",isPET(object),"\n")
    cat("Fragment length:",fragLen(object),"\n")
    cat("Max Bandwidth:", maxBandwidth(object),"\n")
    cat("Using reads files:\n")
    cat(file(object),sep = "\n")
    len = length(seqlengths(regions(object)))
    if( len > 0){
      cat("Using regions for",len,"chromosomes\n")
      show(regions(object))      
    }else{
      cat("**Not regions loaded**\n")
    }
#    cat("---------------------------\n")
})

#' @rdname methods-segvis-readsF
#' @name readsF
setReplaceMethod("readsF",
  signature = signature(object = "segvis",value = "list"),
  definition = function(object,value){
    readsF(object@reads) = value
    return(object)
  }
)                 

#' @rdname methods-segvis-readsR
#' @name readsR
setReplaceMethod("readsR",
  signature = signature(object = "segvis",value = "list"),
  definition = function(object,value){
    readsR(object@reads) = value
    return(object)
  }
)                 



#' @rdname segvis-loadReads
#' @name loadReads
setMethods("loadReads",
  signature = signature(object = "segvis",mc = "numeric"),
  definition = function(object,mc ){

    ## reads the bam file
    message("Reading ",file(object))
    if(isPET(object)){
      message("Setting PET flag")
      pet_flag = scanBamFlag(isPaired = TRUE)
      param = ScanBamParam(which = regions(object),flag = pet_flag)
    }else{
      param = ScanBamParam(which = regions(object))
    }    
    greads = readGAlignmentsFromBam(file(object),
      param = param,use.names = FALSE)
    greads = .data.table.GRanges(as(greads, "GRanges"))
    setkey(greads,seqnames,strand)    
    message("Bam file loaded")

    ## validates the chromosomes, that must coincide between
    ## regions, reads and the ones supplied by user
    chr_reads = unique(greads$seqnames)
    chr_in =sapply(chr_reads,function(x)x%in%chr(object))
    if(!all(chr_in)){
      warning("Removing reads that aren't in ",chr(object))
    }
    greads = greads[chr(object)]

    ## separates by chromosome and strand
    message("Separating reads by strand")
    message("Forward strand reads extracted")
    greads1 = separate.by.chrom(greads,chr(object),"+",mc,sort=TRUE)
    names(greads1) = chr(object)
    message("Reverse strand reads extracted")
    greads2 = separate.by.chrom(greads,chr(object),"-",mc,sort=TRUE)
    names(greads2) = chr(object)
    message("Finished separating reads")

    ## Creates the reads object
    object@reads = new("reads",readsF = greads1,readsR = greads2)
    object@.haveReads = TRUE
    message("Reading bam files... Done")
    return(object)
})

#' @rdname segvis-matchReads
#' @name matchReads
setMethods("matchReads",
  signature = signature(object = "segvis",mc = "numeric"),
  definition = function(object,mc = 8){
    if(object@.haveReads & object@.haveRegions){

      ## separates the regions by chromose, extend the regions by (maxBw - 1)/2 to each
      ## side 
      side = (maxBandwidth(object)-1)/2      
      chr = names(seqlengths(regions(object)))
      match_regions = regions(object)
      start(match_regions) = start(match_regions) - side
      end(match_regions) = end(match_regions) + side      
      match_regions = separate.by.chrom(.data.table.GRanges(match_regions),
        chr, "*",mc,sort=FALSE)
      names(match_regions) = chr

      ## find overlaps between reads and regions   
      overlapsF = mcmapply(.find.overlaps,readsF(object),
        match_regions,SIMPLIFY=FALSE,mc.cores =mc,mc.silent =TRUE)
      overlapsR = mcmapply(.find.overlaps,readsR(object),
        match_regions,SIMPLIFY=FALSE,mc.cores =mc,mc.silent =TRUE)
              
      ## matching reads
      message("Matching reads for forward strand")
      readsF(object) = mcmapply(.match.reads,readsF(object),overlapsF,
        SIMPLIFY= FALSE,mc.cores = mc ,mc.silent =TRUE)
      readsR(object) = mcmapply(.match.reads,readsR(object),overlapsR,
        SIMPLIFY= FALSE,mc.cores = mc ,mc.silent =TRUE)

      object@.readsMatched = TRUE
     return(object)      
    }else{
      warning("Check that both reads and regions are loaded")
    }
})

#' @rdname segvis-getCoverage
#' @name getCoverage
setMethods("getCoverage",
  signature = signature(object = "segvis",mc = "numeric"),
  definition = function(object, mc = 8){
    if(object@.readsMatched == TRUE){

      # init coverage calculation
      chr = names(seqlengths(regions(object)))
      matched_regions = separate.by.chrom(.data.table.GRanges(regions(object)),
          chr, "*",mc,sort=FALSE)
      nregions = mclapply(matched_regions,nrow,mc.cores = mc,mc.silent =TRUE)

      # coverage calculation ,
      message("Calculating coverage")      
      curves = mapply(calculate_chrom_coverage,chr,nregions,
        MoreArgs = list(object,mc),SIMPLIFY=FALSE)      
      names(curves) = chr
      object@profileCurve = curves
      object@.coverageCalculated = TRUE
      message("Coverage done")
      return(object)                
      ## ## ll = lapply(chr,function(chrom,reg){
      ## ##   length(subset(reg,subset = as.character(seqnames(reg)) == chrom))
      ## ## },regions(object))
      ## ## names(ll) = chr
      ## curve = lapply(chr,function(chrom,object,ll,mc){
      ##   message("Retrieving reads for ",chrom)
      ##   l = ll[[chrom]]
      ##   r1 = readsF(object)[[chrom]]
      ##   r2 = readsR(object)[[chrom]]
      ##   m1 = matchF(object)[[chrom]]
      ##   m2 = matchR(object)[[chrom]]
      ##   message("Calculating coverage for ",chrom)
      ##   z = mclapply(1:l,function(i,r1,r2,m1,m2){
      ##     if(isPET(object) | fragLen(object) == 0 ){
      ##       cover = coverage(c(r1[m1[[i]]],r2[m2[[i]]]))[[chrom]]
      ##       }else{
      ##       cover = coverage(c(resize(r1[m1[[i]]],fragLen(object)),
      ##        resize(r2[m2[[i]]],fragLen(object))))[[chrom]]
      ##     }
      ##   return(cover)},r1,r2,m1,m2,mc.cores = mc)        
      ##   message("Coverage calculated for ",chrom)
      ##   return(z)},
      ##   object,ll,mc)     
    }else{
      warning("The reads haven't been matched yet")
    }
})    

# @rdname profile-methods
# @name buildProfileMatrix
# @aliases profile
setMethods("buildProfileMatrix",
  signature = signature(object = "profile",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){
  if(length(unique(width(regions(object))))>1){
      stop("The width of the regions isn't unique, can't build profile matrix")
  }  
  message("Calculating profiles for each chromosome")
  matList = .calculateMAprofile(object,bw,mc)
  message("Joining profiles into matrix")
  matList = lapply(matList,function(x)do.call(rbind,x))
  matList = do.call(rbind,matList)
  message("Profile matrix built")
  return(matList)
})    

# @rdname profile-methods
# @name findSummit
# @aliases profile
setMethods("findSummit",          
  signature = signature(object = "profile",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){
  message("Calculating profiles for each chromosome")
  matList = .calculateMAprofile(object,bw,mc)
  message("Finding summits for each chromosome")
  chr = names(seqlengths(regions(object)))
  if(length(remChr(object)>1))
  {                
    chr = chr[!chr %in% remChr(object)]       
  }else{
    if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
  }
  regions = lapply(chr,function(chrom,reg){
        subset(reg,subset = as.character(seqnames(reg)) == chrom)
      },regions(object))
  names(regions) = chr
  summits_chr = lapply(chr,function(chrom,regions,matList){
    message("Finding summit for regions of ",chrom)    
    reg = regions[[chrom]]
    ll = length(reg)
    rStart = start(reg)
    rEnd = end(reg)
    mat = matList[[chrom]]
    z = mclapply(1:ll,function(i,mat,rStart,rEnd){
      if(all(is.na(mat[[i]]))){
        summit = NA
      }else{
        x = seq(rStart[i],rEnd[i],by = 1)
        summit = x[which.max(mat[[i]])]
      }
      return(summit)
    },mat,rStart,rEnd)
    return(unlist(z))         
  },regions,matList)
  return(unlist(summits_chr))    
})

# @rdname profile-methods
# @name normConst
# @aliases profile
setMethods("countReads",
  signature = signature(object = "profile"),
  definition = function(object){   
    reader = bamReader(file(object),idx=TRUE)
    readMatrix = bamCountAll(reader)
    counts = sum(readMatrix$nAligns)
    return(counts)  
})

# @rdname profile-methods
# @name ProfileMatrix
# @aliases profile
setMethods("ProfileMatrix",
  signature = signature(object = "profile",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc){    
    mat = buildProfileMatrix(object,bw,mc)
    nc = countReads(object)
    nm = name(object)
    gr = regions(object)
    pmatrix = new("profileMatrix",name = nm,regions = gr,profileMat = mat,bandwidth = bw,normConst = nc)
    return(pmatrix)
})    
