
# Methods for profile class

## Get methods

# @rdname profile-methods
# @name name
# @aliases profile
setMethods("name",
  signature = signature(object = "profile"),
  definition = function(object)object@name
)           

# @rdname profile-methods
# @name regions
# @aliases profile
setMethods("regions",
  signature = signature(object = "profile"),
  definition = function(object)object@regions
)           

# @rdname profile-methods
# @name file
# @aliases profile
setMethods("file",
  signature = signature(object ="profile"),
  definition = function(object)object@file
)           

# @rdname profile-methods
# @name fileFormat
# @aliases profile
setMethods("fileFormat",
  signature = signature(object = "profile"),
  definition = function(object)object@fileFormat
)           

# @rdname profile-methods
# @name maxBandwidth
# @aliases profile
setMethods("maxBandwidth",
  signature = signature(object = "profile"),
  definition = function(object)object@maxBandwidth
)           

# @rdname profile-methods
# @name fragLen
# @aliases profile
setMethods("fragLen",
  signature = signature(object = "profile"),
  definition = function(object)object@fragLen
)           

# @rdname profile-methods
# @name remChr
# @aliases profile
setMethods("remChr",
  signature = signature(object = "profile"),
  definition = function(object)object@remChr
)           

# @rdname profile-methods
# @name reads1
# @aliases profile
setMethods("reads1",
  signature = signature(object = "profile"),
  definition = function(object)reads1(object@reads)
)           

# @rdname profile-methods
# @name reads2
# @aliases profile
setMethods("reads2",
  signature = signature(object = "profile"),
  definition = function(object)reads2(object@reads)
)

# @rdname profile-methods
# @name match1
# @aliases profile
setMethods("match1",
  signature = signature(object = "profile"),
  definition = function(object)match1(object@match)
)

# @rdname profile-methods
# @name match2
# @aliases profile
setMethods("match2",
  signature = signature(object = "profile"),
  definition = function(object)match2(object@match)
)

# @rdname profile-methods
# @name profileCurve
# @aliases profile
setMethods("profileCurve",
  signature = signature(object = "profile"),
  definition = function(object)object@profileCurve
)  

## Set methods

# @rdname profile-methods
# @name setName
# @aliases profile
setMethods("setName",
  signature = signature(object = "profile",newName = "character"),
  definition = function(object,newName){
    object@name = newName
    return(object)
})
   
# @rdname profile-methods
# @name setRegions
# @aliases profile
setMethods("setRegions",
  signature = signature(object = "profile",newRegions = "GRanges"),
  definition = function(object,newRegions){
    stopifnot(class(newRegions) == "GRanges")
    object@regions = newRegions
    object@.haveRegions = TRUE
    return(object)
})

# @rdname profile-methods
# @name setMaxBandwidth
# @aliases profile
setMethods("setMaxBandwidth",
  signature = signature(object = "profile", newMaxBandwidth = "numeric"),
  definition = function(object,newMaxBandwidth){
    stopifnot(newMaxBandwidth >= 1)
    stopifnot(newMaxBandwidth %% 2 == 0)
    object@maxBandwidth = newMaxBandwidth
    return(object)
})    

# @rdname profile-methods
# @name setFragLen
# @aliases profile
setMethods("setFragLen",
  signature = signature(object = "profile",newFragLen = "numeric"),
  definition = function(object,newFragLen){
    stopifnot(newFragLen >= 0)
    stopifnot(newFragLen == floor(newFragLen))
    object@fragLen = newFraLen
    return(object)
})    

# @rdname profile-methods
# @name setRemChr
# @aliases profile
setMethods("setRemChr",
  signature = signature(object = "profile",newRemChr = "character"),
  definition = function(object,newRemChr){
    stopifnot(is.character(newRemChr))
    object@remChr = newRemChr
    return(object)
})



# @rdname profile-methods
# @name show
# @aliases profile
setMethods("show",
  signature = signature(object = "profile"),
  definition = function(object){
    cat("---------------------------\n")
    cat("Profile for",name(object),"regions\n")
    cat("Fragment length:",fragLen(object),"\n")
    cat("Max Bandwidth:", maxBandwidth(object),"\n")
    if( length(regions(object)) > 0){
      cat("Using regions for",length(regions(object)),"chromosomes\n")
    }else{
      cat("**Not regions loaded**\n")
    }     
    cat("Using reads files:\n")
    cat(file(object),sep = "\n")
    cat("---------------------------\n")
})
          
# @rdname profile-methods
# @name loadReads
# @aliases profile
setMethods("loadReads",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object,mc ){
    if(fileFormat(object) == "bam"){
      chr = names(seqlengths(regions(object)))     
      if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
      message("Reading ",file(object))
      greads = readGAlignmentsFromBam(file(object),param = NULL,use.names = FALSE)
      greads = as(greads, "GRanges")
      seqlevels(greads) = seqlevelsInUse(greads)
      seqlevels(greads) =seqlevels(regions(object))
      message("Bam file loaded")
      message("Separating by chromosome")      
      greads = mclapply(chr,function(i,greads)
        subset(greads, subset = as.character(seqnames(greads)) == i),greads,mc.cores = mc)
      names(greads) = chr
      message("Separating reads by strand")
      gr1 = mclapply(greads,function(x)sort_by_strand(subset(x,subset = as.character(strand(x))=="+"),"+"),
        mc.cores = mc)
      message("Forward strand reads extracted")
      gr2 = mclapply(greads,function(x)sort_by_strand(subset(x,subset = as.character(strand(x))=="-"),"-"),
        mc.cores = mc)
      message("Reverse strand reads extracted")
      message("Finished separating reads")
      gr1 = GRangesList(gr1)
      gr2 = GRangesList(gr2)
      object@reads = new("reads",reads1 = gr1,reads2 = gr2)
      object@.haveReads = TRUE
      message("Reading bam files... Done")
      return(object)
    }else{# Change warning
      warning("loadReads method only defined for bam file format")
    }
})

# @rdname profile-methods
# @name matchReads
# @aliases profile
setMethods("matchReads",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object,mc = 8){
    if(object@.haveReads & object@.haveRegions){    
      side = (maxBandwidth(object)-1)/2      
      chr = names(seqlengths(regions(object)))
      regions = lapply(chr,function(chrom,reg){
        subset(reg,subset = as.character(seqnames(reg)) == chrom)
      },regions(object))
      regions = GRangesList(lapply(regions,function(x)
        trim(GRanges(seqnames = seqnames(x),ranges = IRanges(start = start(x) - side,
          end = end(x) + side),strand = strand(x)))))
      names(regions) = chr      
      if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]     
      message("Matching reads for forward strand")
      m1 = mclapply(chr,function(chrom,reads,reg,object){
        message("Matching forward reads for ",chrom)
        overlaps = findOverlaps(reg[[chrom]],resize(reads[[chrom]],fragLen(object)))       
        mm = lapply(1:length(reg[[chrom]]),
          function(i,overlaps)subjectHits(subset(overlaps,subset = queryHits(overlaps) == i)),
          overlaps)  
        message("Forward strand matching for ",chrom," done");return(mm)},
        reads1(object),regions,object,mc.cores = mc)
      names(m1) = chr            
      message("Forward strand done")
      message("Matching reads for reverse strand")
      m2 = mclapply(chr,function(chrom,reads,reg,object){
        message("Matching reverse reads for ",chrom)
        overlaps = findOverlaps(reg[[chrom]],resize(reads[[chrom]],fragLen(object)))
        mm = lapply(1:length(reg[[chrom]]),
          function(i,overlaps)subjectHits(subset(overlaps,subset = queryHits(overlaps) == i)),
          overlaps)
        message("Reverse strand matching for ",chrom," done");return(mm)},
        reads2(object),regions,object,mc.cores = mc)
      names(m2) = chr
      message("Reverse strand done")
      object@match = new("match",match1 = m1,match2 = m2)
      object@.readsMatched = TRUE
     return(object)      
    }else{
      warning("Check that both reads and regions are loaded")
    }
})

# @rdname profile-methods
# @name getCoverage
# @aliases profile
setMethods("getCoverage",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object, mc = 8){
    if(object@.readsMatched == TRUE){      
      chr = names(seqlengths(regions(object)))
      if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]     
      ll = lapply(chr,function(chrom,reg){
        length(subset(reg,subset = as.character(seqnames(reg)) == chrom))
      },regions(object))
      names(ll) = chr           
      curve = lapply(chr,function(chrom,object,ll,mc){
        message("Retrieving reads for ",chrom)        
        l = ll[[chrom]]
        r1 = reads1(object)[[chrom]]
        r2 = reads2(object)[[chrom]]
        m1 = match1(object)[[chrom]]
        m2 = match2(object)[[chrom]]
        message("Calculating coverage for ",chrom)
        z = mclapply(1:l,function(i,r1,r2,m1,m2){
          coverage(c(resize(r1[m1[[i]]],fragLen(object)),
                     resize(r2[m2[[i]]],fragLen(object))))[[1]]}
                 ,r1,r2,m1,m2,mc.cores = mc)        
        message("Coverage calculated for ",chrom)
        return(z)},
        object,ll,mc)     
      names(curve) = chr
      object@profileCurve = curve
      object@.coverageCalculated = TRUE
      message("Coverage done")
      return(object)     
    }else{
      warning("The reads haven't been matched yet")
    }
})    

# @rdname profile-methods
# @name buildProfileMat
# @aliases profile
setMethods("buildProfileMatrix",
  signature = signature(object = "profile",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){
  if(object@.coverageCalculated){    
    ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)
    chr = names(seqlengths(regions(object)))
    if(remChr(object) != ""){
      chr = chr[!chr %in% remChr(object)]
    }
    if(length(unique(width(regions(object))))>1){
      stop("The width of the regions isn't unique, can't build profile matrix")}
    side = (maxBandwidth(object)-1)/2
    regions = GRangesList(lapply(chr,function(x,regions)
      subset(regions,subset = as.character(seqnames(regions)) == x),regions(object)))
    names(regions) =chr
    matList = lapply(chr,function(chrom,object,regions,side,mc){
      message("Calculating profile for ",chrom)
      ll = length(regions[[chrom]])
      regionStart = start(regions[[chrom]])-side
      regionEnd = end(regions[[chrom]])+side
      stepList = profileCurve(object)[[chrom]]   
      mclapply(1:ll,function(i,regionStart,regionEnd,stepList,bw,side){        
        z = stepList[[i]]
        x = seq(regionStart[i],regionEnd[i],by=1)
        if(nrun(z)==1){
          y = rep(runValue(z),length(x))
        }else{
          xp = cumsum(runLength(z)[1:(nrun(z)-1)])       
          yp = runValue(z)        
          y = stepfun(xp,yp)(x)
          y =  ma(y,bw)
        }              
        y = y[-c(1:side)]
        y = y[1:(length(y) - side)]                            
    return(y)},regionStart,regionEnd,stepList,bw,side,mc.cores = mc)
      },object,regions,side,mc)
    matList = lapply(matList,function(x)do.call(rbind,x))
    matList = do.call(rbind,matList)
    message("Profile matrix built")
    return(matList)
  }else{
    warning("The coverage haven't been calculated yet")
  }
})    
