#' @import ggplot2
NULL

#' @name plot_profiles
#' @rdname plot_profiles
#' @aliases plot_profiles,ANY-method
setMethods("plot_profiles",
  signature = signature(object = "segvis_block_list",
    condition = "ANY",coord = "numeric",FUN = "function",...="ANY",
    mc = "numeric"),
  definition = function(object,condition,coord,FUN,...,mc){

    if(is.null(names(object))){nms = 1:length(object)}else{nms = names(object)}
    if(!missing(condition)){
      conds <- mclapply(object,.subset_logical,substitute(condition),mc.cores = mc)
    }else{      
      nregions <- lapply(object,function(x)length(regions(x)))
      conds <- lapply(nregions,function(x)rep(TRUE,x))
    }    
    subsets <- mcmapply(.filter_sb,object,conds,SIMPLIFY=FALSE,mc.cores = mc)
    widths <- mclapply(subsets,function(x)width(regions(x)),mc.cores = mc)
    if(length(unique(unlist(widths)))>1){
      stop("The supplied regions doesn't have the same width")
    }else      
    plot_width <- unique(unlist(widths))
    if(missing(coord)){
      coord <- 1:plot_width
    }    
    profiles <- mclapply(subsets,summarize,FUN,...,mc.cores = mc)
    plot_data <- mcmapply(create_plot_data,profiles,nms,MoreArgs = list(coord),
      SIMPLIFY = FALSE,mc.silent=TRUE,mc.cores = mc,mc.preschedule=TRUE)
    plot_data <- do.call(rbind,plot_data)
    x <- y <- NULL
    out <- ggplot(plot_data,aes(x,y,colour = as.factor(condition)))+geom_line(size=1.1)+
      scale_colour_discrete(guide = guide_legend(title = "condition"))+
      xlab("genomic coordinates")+ylab("normalized coverage")        
    return(out)                  
})

#' @name plot_data
#' @rdname plot_data
#' @aliases plot_data,ANY-method
setMethods("plot_data",
  signature = signature(object = "segvis_block_list",
    condition = "ANY",coord = "numeric",FUN = "function",...="ANY",
    mc = "numeric"),
  definition = function(object,condition,coord,FUN,...,mc){

  if(is.null(names(object))){nms = 1:length(object)}else{nms = names(object)}
  if(!missing(condition)){
    conds <- mclapply(object,.subset_logical,substitute(condition),mc.cores = mc)
  }else{      
    nregions <- lapply(object,function(x)length(regions(x)))
    conds <- lapply(nregions,function(x)rep(TRUE,x))
  }    
  subsets <- mcmapply(.filter_sb,object,conds,SIMPLIFY=FALSE,mc.cores = mc)
  widths <- mclapply(subsets,function(x)width(regions(x)),mc.cores = mc)
  if(length(unique(unlist(widths)))>1){
    stop("The supplied regions doesn't have the same width")
  }else      
  plot_width <- unique(unlist(widths))
  if(missing(coord)){
    coord <- 1:plot_width
  }
  profiles <- mclapply(subsets,summarize,FUN,...,mc.cores = mc)
  plot_data <- mcmapply(create_plot_data,profiles,nms,MoreArgs = list(coord),
    SIMPLIFY = FALSE,mc.silent = TRUE,mc.cores = mc,mc.preschedule=TRUE)
  plot_data <- do.call(rbind,plot_data)    
 return(plot_data)
})

