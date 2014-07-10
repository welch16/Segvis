
# Methods for profile class

## Get methods

#' @rdname profile-methods
#' @name name
#' @aliases profile
setMethods("name",
  signature = signature(object = "profile"),
  definition = function(object)object@name
)           

#' @rdname profile-methods
#' @name regions
#' @aliases profile
setMethods("regions",
  signature = signature(object = "profile"),
  definition = function(object)object@regions
)           

#' @rdname profile-methods
#' @name files
#' @aliases profile
setMethods("files",
  signature = signature(object ="profile"),
  definition = function(object)object@files
)           

#' @rdname profile-methods
#' @name fileFormat
#' @aliases profile
setMethods("fileFormat",
  signature = signature(object = "profile"),
  definition = function(object)object@fileFormat
)           

#' @rdname profile-methods
#' @name maxBandwidth
#' @aliases profile
setMethods("maxBandwidth",
  signature = signature(object = "profile"),
  definition = function(object)object@maxBandwidth
)           

#' @rdname profile-methods
#' @name fragLen
#' @aliases profile
setMethods("fragLen",
  signature = signature(object = "profile"),
  definition = function(object)object@fragLen
)           

#' @rdname profile-methods
#' @name remChr
#' @aliases profile
setMethods("remChr",
  signature = signature(object = "profile"),
  definition = function(object)object@remChr
)           

#' @rdname profile-methods
#' @name readsList
#' @aliases profile
setMethods("readsList",
  signature = signature(object = "profile"),
  definition = function(object)object@readsList           
)
           
#' @rdname profile-methods
#' @name matchList
#' @aliases profile
setMethods("matchList",
  signature = signature(object = "profile"),
  definition = function(object)object@matchList
)

#' @rdname profile-methods
#' @name profileCurve
#' @aliases profile
setMethods("profileCurve",
  signature = signature(object = "profile"),
  definition = function(object)object@profileCurve
)  

## Set methods

#' @rdname profile-methods
#' @name setName
#' @aliases profile
setMethods("setName",
  signature = signature(object = "profile",newName = "character"),
  definition = function(object,newName){
    object@name = newName
    return(object)
})

#' @rdname profile-methods
#' @name setRegions
#' @aliases profile
setMethods("setRegions",
  signature = signature(object = "profile",newRegions = "GRangesList"),
  definition = function(object,newRegions){
    stopifnot(class(newRegions) == "GRangesList")
    object@regions = newRegions
    object@.haveRegions = TRUE
    return(object)
})

#' @rdname profile-methods
#' @name setMaxBandwidth
#' @aliases profile
setMethods("setMaxBandwidth",
  signature = signature(object = "profile", newMaxBandwidth = "numeric"),
  definition = function(object,newMaxBandwidth){
    stopifnot(newMaxBandwidth >= 1)
    stopifnot(newMaxBandwidth %% 2 == 0)
    object@maxBandwidth = newMaxBandwidth
    return(object)
})    

#' @rdname profile-methods
#' @name setFragLen
#' @aliases profile
setMethods("setFragLen",
  signature = signature(object = "profile",newFragLen = "numeric"),
  definition = function(object,newFragLen){
    stopifnot(newFragLen >= 0)
    stopifnot(newFragLen == floor(newFragLen))
    object@fragLen = newFraLen
    return(object)
})    

#' @rdname profile-methods
#' @name setRemChr
#' @aliases profile
setMethods("setRemChr",
  signature = signature(object = "profile",newRemChr = "character"),
  definition = function(object,newRemChr){
    stopifnot(is.character(newRemChr))
    object@remChr = newRemChr
    return(object)
})
              
#' @rdname profile-methods
#' @name show
#' @aliases profile
setMethods("show",
  signature = signature(object = "profile"),
  definition = function(object){
    cat("---------------------------\n")
    cat("Profile for",name(object),"peaks\n")
    cat("Fragment length:",fragLen(object),"\n")
    cat("Max Bandwidth:", maxBandwidth(object),"\n")
    if(q <- length(regions(object)) > 0){
      cat("Using regions for",length(regions(object)),"chromosomes\n")
    }else{
      cat("**Not regions loaded**\n")
    }     
    cat("Using reads files:\n")
    cat(files(object),sep = "\n")
    cat("---------------------------\n")
})
          
#' @rdname profile-methods
#' @name loadReads
#' @aliases profile
setMethods("loadReads",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object,mc = 8){
    if(fileFormat(object) == "bam"){
      chr = names(seqlengths(regions(object)))
      if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
      message("Starting to read bam files")
      greads = lapply(files(object),FUN = readGAlignmentsFromBam,param = NULL,use.names = FALSE)      
      greads = lapply(greads,FUN = as, "GRanges")
      message("Bam files loaded")
      message("Separating by chromosome")
      greads = lapply(greads,function(x,chr,mc){
        z =mclapply(chr,function(i,x)subset(x, subset = as.character(seqnames(x)) == i),x,mc.cores = mc);
        names(z) = chr;
        return(z)
        },chr,mc)
      message("Separating reads by strand")
      gr1 = lapply(greads,function(x,mc)
          mclapply(x,function(y)sort_by_strand(subset(y,subset = as.character(strand(y)) == "+"),"+"),
                 mc.cores = mc),mc)      
      gr2 = lapply(greads,function(x,mc)
          mclapply(x,function(y)sort_by_strand(subset(y,subset = as.character(strand(y)) == "-"),"-"),
                mc.cores = mc),mc)     
      gr1 = lapply(gr1,FUN = GRangesList)
      gr2 = lapply(gr2,FUN = GRangesList)
      object@readsList = lapply(1:length(gr1),function(i,gr1,gr2)
        new("reads",reads1 = gr1[[i]],reads2 = gr2[[i]]),gr1,gr2)
      object@.haveReads = TRUE
      message("Reading bam files... Done")
      return(object)
    }else{
      warning("loadReads method only defined for bam file format")
    }
})

#' @rdname profile-methods
#' @name matchReads
#' @aliases profile
setMethods("matchReads",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object,mc = 8){
    if(object@.haveReads & object@.haveRegions){
      side = (maxBandwidth(object)-1)/2
      chr = names(seqlengths(regions(object)))
      if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
      message("Matching reads for + strand")
      m1 = lapply(readsList(object),function(x,chr,object,mc){
        z = mclapply(chr,function(chrom,x,object){          
          message(chrom," started");          
          mm1 =match_reads(start(regions(object)[[chrom]])-side,
          end(regions(object)[[chrom]])+side,start(reads1(x)[[chrom]]),end(reads1(x)[[chrom]]),
          as.character(strand(reads1(x)[[chrom]])),fragLen(object));
          message(chrom," done");return(mm1)},x,object, mc.cores=mc)
        names(z) = chr
        return(z)},chr,object,mc)
      message("+ strand done")
      message("Matching reads for - strand")
      m2 = lapply(readsList(object),function(x,chr,object,mc){        
        z = mclapply(chr,function(chrom,x,object){message(chrom," started");
          mm2=match_reads(start(regions(object)[[chrom]])-side,end(regions(object)[[chrom]])+side,
            start(reads2(x)[[chrom]]),end(reads2(x)[[chrom]]),
          as.character(strand(reads2(x)[[chrom]])),fragLen(object));
          message(chrom," done");return(mm2)},x,object, mc.cores=mc)
        names(z) = chr
        return(z)},chr,object,mc)
      message("- strand done")
      object@matchList = lapply(1:length(m1),function(i,m1,m2)
        new("match",match1 = m1[[i]],match2 = m2[[i]]),m1,m2)
      object@.readsMatched = TRUE
     return(object)      
    }else{
      warning("Check that both reads and regions are loaded")
    }
})

#' @rdname profile-methods
#' @name getCoverage
#' @aliases profile
setMethods("getCoverage",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object, mc = 8){
    if(object@.readsMatched == TRUE){
      chr = names(seqlengths(regions(object)))
      if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
      object@profileCurve = lapply(chr,function(chrom,object,mc){message("Retrieving reads for ",chrom)
         ll = length(regions(object)[[chrom]])
         r1 = reads1(readsList(object)[[1]])[[chrom]]
         r2 = reads2(readsList(object)[[1]])[[chrom]]
         mclapply(1:ll,function(i,object,chrom,r1,r2){
         z = coverage(c(r1[ match1(matchList(object)[[1]])[[chrom]] [[i]] ],
                        r2[ match2(matchList(object)[[1]])[[chrom]] [[i]] ]))[[chrom]]
      return(z)},object,chrom,r1,r2,mc.cores = mc)},object,mc)
      names(object@profileCurve) = chr
      object@.coverageCalculated = TRUE
      message("Coverage done")
      return(object)     
    }else{
      warning("The reads haven't been matched yet")
    }
})    

#' @rdname profile-methods
#' @name buildProfileMat
#' @aliases profile
setMethods("buildProfileMat",
  signature = signature(object = "profile",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){
  if(object@.coverageCalculated){
    ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)
    chr = names(seqlengths(regions(object)))
    if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
    side = (maxBandwidth(object)-1)/2
    matList = lapply(chr,function(chrom,object,side,mc){
      message("Calculating profile for ",chrom)
      ll = length(regions(object)[[chrom]])
      regionStart = start(regions(object)[[chrom]])-side
      regionEnd = end(regions(object)[[chrom]])+side
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
      },object,side,mc)
    matList = lapply(matList,function(x)do.call(rbind,x))
    matList = do.call(rbind,matList)
    message("Profile matrix built")
    return(matList)
  }else{
    warning("The coverage haven't been calculated yet")
  }
})    
