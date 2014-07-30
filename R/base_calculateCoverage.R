# Calculates the moving average of a veces
# @param x - Numeric vector
# @param n - Integer, represents the bandwidth of the moving average
# @return A Numeric Vector
.ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)

# Calculates the MA average profile without width restriction
# The idea of this funcion is to use it to plot peaks, find summit or build the profile matrix
# @param object - Profile object
# @param bw - Numeric bandwidth to be used
# @param mc - Numeric cores to be used
# @return List - Is a list for each chromosome with an element for each region
.calculateMAprofile <- function(object,bw,mc)
{  
  if(object@.coverageCalculated){
  stopifnot(bw <= maxBandwidth(object))
  chr = names(seqlengths(regions(object)))
  if(length(remChr(object)>1))
  {                
    chr = chr[!chr %in% remChr(object)]       
  }else{
    if(remChr(object) != "")chr = chr[!chr %in% remChr(object)]
  }
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
    mclapply(1:ll,function(i,regionStart,regionEnd,stepList,bw,side){       z = stepList[[i]]
      x = seq(regionStart[i],regionEnd[i],by=1)
      if(nrun(z)==1){
        if(runValue(z) == 0){
          y = rep(NA,length(x))
        }else{
          y = rep(runValue(z),length(x))
        }
      }else{
        xp = cumsum(runLength(z)[1:(nrun(z)-1)])       
        yp = runValue(z)        
        y = stepfun(xp,yp,right = TRUE)(x)
        y =  .ma(y,bw)
      }              
      y = y[-c(1:side)]
      y = y[1:(length(y) - side)]                            
      return(y)},regionStart,regionEnd,stepList,bw,side,mc.cores = mc)
    },object,regions,side,mc)
  names(matList) = chr
  return(matList)
  }else{
    warning("The coverage haven't been calculated yet")
  }
}    
