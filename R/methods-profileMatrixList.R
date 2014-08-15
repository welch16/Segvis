 
# @rdname profileMatrixList-methods
# @name plot.profiles
setMethods("plot.profiles",
  signature = signature(object = "profileMatrixList",condition = "ANY",coord = "numeric",trim = "numeric"),
  definition = function(object,condition,coord,trim){    
    ncols = unique(sapply(object,function(x)ncol(profileMat(x))))    
    if(length(ncols)>1)stop("The profile matrices doesn't have the same number of columns")
    if(missing(coord))coord = 1:ncols    
    if(missing(trim))trim = 0    
    stopifnot(length(coord) == ncols)   
    if(!missing(condition)){
      ll = list()
      if(is.null(names(object))){
        names(object) = 1:length(object)
      }    
      for(i in names(object)){        
        x = object[[i]]
        cond = .subset_profileMat_logical(x,substitute(condition))
        ll[[i]] = .filter_profileMat(x,cond)
      }     
    }else{
      ll = as(object, "list")
    }    
    pairs = lapply(ll,function(x,coord)data.frame(coord = coord,profile = meanProfile(x)),coord)
    pairs = lapply(names(pairs),function(x,pairs){
      pairs[[x]]$group = x
      return(pairs[[x]])},pairs)
    data = do.call(rbind,pairs)  
    data$group = factor(data$group,levels = names(object))
    out = ggplot(data,aes(coord,profile,colour = group,linetype = group))+geom_line(size = 1)
    return(out)
})    


