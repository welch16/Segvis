
# Methods for profile class

## Get methods

#' name
#'
#' @param profile object
#' @return character. The name of the object
#' @docType methods
#' @rdname profile-methods
setMethods("name",
  signature = signature(object = "profile"),
  definition = function(object)object@name
)           

#' regions
#'
#' @param profile object
#' @return GRangesList. The regions used to calculate the coverage plots, as a GRangesList separated by chromosome
#' @docType methods
#' @rdname profile-methods
setMethods("regions",
  signature = signature(object = "profile"),
  definition = function(object)object@regions
)           

#' files
#'
#' @param profile object
#' @return chracter. A chracter vector with the names of the bedfiles used to create the readsList object
#' @docType methods
#' @rdname profile-methods
setMethods("files",
  signature = signature(object ="profile"),
  definition = function(object)object@files
)           

#' fileFormat
#'
#' @param profile object
#' @return character. A character with the file format used
#' @docType methods
#' @rdname profile-methods
setMethods("fileFormat",
  signature = signature(object = "profile"),
  definition = function(object)object@fileFormat
)           

#' maxBandwidth
#'
#' @param profile object
#' @return numeric. A number with the max bandwidth possible to smooth the profiles
#' @docType methods
#' @rdname profile-methods
setMethods("maxBandwidth",
  signature = signature(object = "profile"),
  definition = function(object)object@maxBandwidth
)           

#' fragLen
#'
#' @param profile object
#' @return numeric. A numeric value representing the width of the extended fragment reads
#' @docType methods
#' @rdname profile-methods
setMethods("fragLen",
  signature = signature(object = "profile"),
  definition = function(object)object@fragLen
)           

#' readsList
#'
#' @param profile object
#' @return list. A list made off reads objects. One list for each replicate. Must coincide with size of bedfiles
#' docType methods
#' @rdname profile-methods
setMethods("readsList",
  signature = signature(object = "profile"),
  definition = function(object)object@readsList           
)
           
#' matchList
#'
#' @param profile object
#' @return list. A list made off match objects. One list for each replicate. Must coincide with size of bedfiles
#' @docType methods
#' @rdname profile-methods
setMethods("matchList",
  signature = signature(object = "profile"),
  definition = function(object)object@matchList
)

#' profileCurve
#'
#' @param profile object
#' @return RleList. A list made of an Rle object for each region
#' @docType methods
#' @rdname profile-methods
setMethods("profileCurve",
  signature = signature(object = "profile"),
  definition = function(object)object@profileCurve
)  

## Set methods

#' setName
#'
#' @param profile object
#' @param newName character
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setName",
  signature = signature(object = "profile",newName = "character"),
  definition = function(object,newName){
    object@name = newName
    return(object)
})

#' setRegions
#'
#' @param profile object
#' @param newRegions GRangesList
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setRegions",
  signature = signature(object = "profile",newRegions = "GRangesList"),
  definition = function(object,newRegions){
    stopifnot(class(newRegions) == "GRangesList")
    object@regions = newRegions
    object@.haveRegions = TRUE
    return(object)
})

#' setMaxBandwidth
#'
#' @param profile object
#' @param newMaxBandwidth Numeric value, must be odd and greater or equal than one
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setMaxBandwidth",
  signature = signature(object = "profile", newMaxBandwidth = "numeric"),
  definition = function(object,newMaxBandwidth){
    stopifnot(newMaxBandwidth >= 1)
    stopifnot(newMaxBandwidth %% 2 == 0)
    object@maxBandwidth = newMaxBandwidth
    return(object)
})    

#' setFragLen
#'
#' @param profile object
#' @param newFragLen Numeric value, must be greater or equal to zero
#' @return profile object
#' @docType methods
#' @rdname profile-methods
setMethods("setFragLen",
  signature = signature(object = "profile",newFragLen = "numeric"),
  definition = function(object,newFragLen){
    stopifnot(newFragLen >= 0)
    stopifnot(newFragLen == floor(newFragLen))
    object@fragLen = newFraLen
    return(object)
})    

#' show
#'
#' @param profile object
#' @docType methods
#' @rdname profile-methods
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
          
#' loadReads
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @docType methods
#' @rdname profile-methods
setMethods("loadReads",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object,mc = 8){
    if(fileFormat(object) == "bam"){
      chr = names(seqlengths(regions(object)))
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

#' matchReads
#' 
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @docType methods
#' @rdname profile-methods
setMethods("matchReads",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object,mc = 8){
    if(object@.haveReads & object@.haveRegions){
      side = (maxBandwidth(object)-1)/2
      chr = names(seqlengths(regions(object)))
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

#' getCoverage
#'
#' @param profile object
#' @param mc, the number of cores used with parallel
#' @docType methods
#' @rdname profile-methods
setMethods("getCoverage",
  signature = signature(object = "profile",mc = "numeric"),
  definition = function(object, mc = 8){
    if(object@.readsMatched == TRUE){
      chr = names(seqlengths(regions(object)))     
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
      return(object)     
    }else{
      warning("The reads haven't been matched yet")
    }
})    

#' buildProfileMat
#'
#' @param profile object
#' @param mc, the number of cores used with parallel
#' @param bw, the bandwidth used to smooth the profile
#' @docType methods
#' rdname methods-profile
setMethods("buildProfileMat",
  signature = signature(object = "profile",bw = "numeric",mc = "numeric"),
  definition = function(object,bw,mc=8){
  if(object@.coverageCalculated){
    ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)
    chr = names(seqlengths(regions(object)))
    side = (maxBandwidth(object)-1)/2
    matList = lapply(chr,function(chrom,object,side,mc){
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
    return(matList)
  }else{
    warning("The coverage haven't been calculated yet")
  }
})    
