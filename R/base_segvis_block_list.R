##' @import ggplot2
NULL

create_plot_data <- function(counts,name,coord)
{
  dt <- data.table(x=coord,y =counts,condition = name)
  return(dt)
}



##' @rdname plot_profiles-methods
##' @name plot_profiles-methods
##' @aliases plot_profiles,ANY-method
##' @docType methods
##' @exportMethod plot_profiles
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


##' Generates the data to be plotted
##'
##' This function generates a data.table with the profiles that would be in the profile plot
##'
##' @param object segvis_block_list object
##'
##' @param condition This is an expression used to filter the \code{segvis_block} objects, if there is no condition, then it doesn't filter the objects
##'
##' @param coord Numeric vector representing the coordinates in the a-axis of the plot, by default uses the natural index \code{1:ncol(segvis_block)}
##'
##' @param mc Numeber of cores used by parallel
##'
##' @param FUN function used to summarize the profiles
##'
##' @param ... Additional arguments of \code{FUN}
##'
##' @export
##'
##' @return A data.table object
##'
##' @rdname plot_data
##' @name plot_data
##' 
##' @seealso \code{\link{segvis_block_list-class}},\code{\link{plot_profiles}}
##'
##' 
##' @examples
##' \dontrun{
##'
##'  plot_data(segvis_block_list,mc = 8, mean, trim = .05)
##'
##'  ## if the segvis_block_list have a column called row
##'  plot_data(segvis_block_list,condition = row == 1,mc = 8,FUN = function(x)x)
##'
##' }
plot_data <- function(object,condition,coord = NULL,mc,FUN,...)
{
  ## stopifnot(class(object) == "segvis_block_list")
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
  if(is.null(coord)){
    coord <- 1:plot_width
  }
  profiles <- mclapply(subsets,summarize,FUN,...,mc.cores = mc)
  plot_data <- mcmapply(create_plot_data,profiles,nms,MoreArgs = list(coord),
    SIMPLIFY = FALSE,mc.silent = TRUE,mc.cores = mc,mc.preschedule=TRUE)
  plot_data <- do.call(rbind,plot_data)    
 return(plot_data)
}

