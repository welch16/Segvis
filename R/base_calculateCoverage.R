# Calculates the moving average of a veces
# @param x - Numeric vector
# @param n - Integer, represents the bandwidth of the moving average
# @return A Numeric Vector
.ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)



.calc_profile <- function(regionStart,regionEnd,step_fn,bw,side)
{
  z = step_fn
  x = seq(regionStart,regionEnd,by=1)
  if(!is.null(z)){
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
    if(side > 0){
      y = y[-c(1:side)]
      y = y[1:(length(y) - side)]
    }
  }else{
    y = rep(0,length(x))
    if(side > 0){
      y = y[-c(1:side)]
      y = y[1:(length(y) - side)]
    }       
  }
  return(as.numeric(y))
}

.find_summit <- function(curve,regionStart,regionEnd)
{
  if(all(is.na(curve))){
    summit = NA
  }else{
    x = seq(regionStart,regionEnd,by=1)
    summit = x[which.max(curve)]
  }
  return(summit)
}

# Calculates the MA average profile without width restriction
# The idea of this funcion is to use it to plot peaks, find summit or build the profile matrix
# @param object - Profile object
# @param bw - Numeric bandwidth to be used
# @param mc - Numeric cores to be used
# @return List - Is a list for each chromosome with an element for each region
.calculate_profile_curves<- function(object,bw,mc)
{  
  if(object@.coverageCalculated){

    stopifnot(bw <= maxBandwidth(object))

    ## init algorithm to calculate curves
    chr = names(seqlengths(regions(object)))
    side = (maxBandwidth(object)-1)/2
    match_regions = separate.by.chrom(.data.table.GRanges(regions(object)),
      chr, "*",mc,sort=FALSE)
    names(match_regions) = chr

    ## calculate profile curves
    curves = lapply(chr,function(chrom,object,regions,side,mc){
      message("Calculating profile for ",chrom)
      chr_reg = regions[[chrom]]
      nreg = length(chr_reg)
      regionStart = chr_reg[,(start)]-side 
      regionEnd = chr_reg[,(end)] + side 
      stepList = profiles(object)[[chrom]]
      chr_curves = mcmapply(.calc_profile,regionStart,regionEnd,stepList,
        MoreArgs = list(bw,side),SIMPLIFY=FALSE,mc.cores=mc,mc.silent=TRUE,
        mc.preschedule =FALSE)
      return(chr_curves)
    },object,match_regions,side,mc)
    names(curves) = chr

    return(curves)
  }else{
    warning("The coverage haven't been calculated yet")
  }
}    


