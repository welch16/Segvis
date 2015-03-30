#' @title Set of functions to subset segvis_block objects
#' @description \code{.subset_logical} Function to obtain a logical vector of the elements in regions(object) such that satisfy condition_all
#' @param object segvis_block object
#' @param condition_call name object with the filtering condition
#' @return A logical vector with the same length as the number of regions
#' @rdname filter_sb
.subset_logical<- function(object,condition_call){
  cond = as.logical(eval(condition_call,as(regions(object),"data.frame"),parent.frame()))
  return(cond)
}

#' @description \code{.filter_sb} Filters the segvis_block object using a logical vector as an input, the cond vector can be obtained by using \code{.subset_logical}
#' @param cond Logical Vector with the peaks that satisfy certain condition
#' @return A segvis_block object with the peaks that such that cond is TRUE
#' @rdname filter_sb
.filter_sb <- function(object,cond){
  if(any(is.na(cond))){
    warning("There are regions impossible to evaluate")
    cond[is.na(cond)] = FALSE
  }
  coverage = copy(cover_table(object))
  lengths = coverage[,length(chr),by=match][,(V1)]
  extended_cond = unlist(mapply(function(x,l)rep(x,l),cond,lengths,SIMPLIFY=FALSE))

  # Add condition to data.table, filter and remove condition
  out_regions = regions(object)[cond]
  rm(cond)
  coverage[,cond:=extended_cond]
  coverage = coverage[cond == TRUE]
  coverage[,cond:=NULL]
    
  out = new("segvis_block",name = name(object),regions = out_regions,
    bandwidth = bandwidth(object),normConst =  normConst(object),
    cover_table = coverage,.isScaled = object@.isScaled)
  
  return(out)    
}
