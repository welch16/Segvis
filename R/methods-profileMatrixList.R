
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
      for(i in names(object)){
        x = object[[i]]
        env = list2env(as(regions(x),"data.frame"),envir = environment())
        cond = eval(substitute(condition),envir = env)
        ll[[i]] = new("profileMatrix",name = name(x),regions = regions(x)[cond],bandwidth = bandwidth(x),normConst = normConst(x),profileMat = profileMat(x)[cond,],.isScaled = x@.isScaled)
      }     
    }else{
      ll = as(object, "list")
    }
    pairs = lapply(ll,function(x,coord)data.frame(coord = coord,profile = meanProfile(x)),coord)
    pairs = lapply(names(pairs),function(x,pairs){
      pairs[[x]]$group = x
      return(pairs[[x]])},pairs)
    data = do.call(rbind,pairs)
    data$group = factor(data$group)
    out = ggplot(data,aes(coord,profile,colour = group))+geom_line(size = 1.2)
    return(out)
})    


