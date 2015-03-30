
# generic methods for various classes

#' @title Get and set methods for segvis class
#'
#' @description This methods determine how to get and set attributes dor the segvis class. The structure is similar to the one use by bioconductor objects
#'
#' @param object An object of class segvis
#'
#' @return Returns the value with the same name as the function
#'
#' @export
#'
#' @docType methods
#'
#' @seealso \code{\link{segvis-class}}
#'
#' @rdname methods-segvis-gs
setGeneric("name",
  function(object)
  standardGeneric("name")
)

#' @param value The class of this parameter depends on the function, for \code{file<-}, \code{name<-} and \code{chr<-} is a string; for \code{fragLen<-} and  \code{maxBandwidth} is a numeric value; for \code{regions<-} is a GRanges object and for \code{isPET<-} is a logical value
#'
#' @return Returns a segvis object with the value indicated by the functions being updated with value
#'
#' @rdname methods-segvis-gs
setGeneric("name<-",
  function(object,value)
  standardGeneric("name<-")
)           

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("regions",
  function(object)
  standardGeneric("regions")
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("regions<-",
  function(object,value)
  standardGeneric("regions<-")
)

# Generic methods for profile class

##  Get methods

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("file",           
  function(object)
  standardGeneric("file")           
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("file<-",
  function(object,value)
  standardGeneric("file<-")
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("maxBandwidth",
  function(object)
  standardGeneric("maxBandwidth")           
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("maxBandwidth<-",
  function(object,value)
  standardGeneric("maxBandwidth<-")
)           

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("fragLen",
  function(object)
  standardGeneric("fragLen")           
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("fragLen<-",
  function(object,value)
  standardGeneric("fragLen<-")
)           

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("chr",
  function(object)
  standardGeneric("chr")
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("chr<-",
  function(object,value)
  standardGeneric("chr<-")           
)

#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("isPET",
  function(object)
  standardGeneric("isPET")           
)
           
#' @export
#' @docType methods
#' @rdname methods-segvis-gs
setGeneric("isPET<-",
  function(object,value)
  standardGeneric("isPET<-")           
)           

#' @title profiles methods for segvis class
#' @description Returns the coverage for each region considered
#' @param object A segvis object
#' @return RleList. A list made of an Rle object for each region
#' @export
#' @docType methods
#' @seealso \code{\link{segvis-class}} and \code{\link{segvis-getCoverage}}
#' @rdname methods-segvis-gs
setGeneric("profiles",
  function(object)
  standardGeneric("profiles")
)  
    
#' @title loadReads method for segvis class
#' @description Load the fragment stored in the file slot of the segvis object. The reads are divided by chromosome and by strand
#' @param object segvis object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{reads-class}}
#' @rdname segvis-loadReads
#' @examples
#' \dontrun{ loadReads(profile_object,mc=8)}
setGeneric("loadReads",
  function(object,mc)
  standardGeneric("loadReads")
)           

#' @title matchReads methods for segvis class
#' @description Match the reads to the extended regions stored in the object slots. 
#' @param object segvis object 
#' @param mc numeric, the number of cores used with parallel
#' @return segvis object
#' @export
#' @docType methods
#' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{reads-class}}
#' @rdname segvis-matchReads
#' @examples
#' \dontrun{ matchReads(segvis_object,mc=8)}
setGeneric("matchReads",
  function(object,mc)
  standardGeneric("matchReads")
)           

#' @title getCoverage method for profile class
#' @description Calculate  the coverage using the reads matched by the matchReads method. It returns a Rle object for each region.
#' @param object segvis object
#' @param mc numeric, the number of cores used with parallel
#' @return segvis object
#' @export
#' @docType methods
#' @seealso \code{\link{loadReads}} and \code{\link{matchReads}}
#' @rdname segvis-getCoverage
#' @examples
#' \dontrun{ getCoverage(segvis_object,mc=8)}
setGeneric("getCoverage",
  function(object,mc)
  standardGeneric("getCoverage")           
)

#' @title joinProfiles method for segvis class
#' @description Joins all the coverages into a data.table with columns chr|match|coord|tagCounts
#' @param object segvis object
#' @param bw Numeric values with the bandwidth used to smooth the profiles, must be and odd number less or equal than the maxBandwidth
#' @param mc Numeric value representing the number of cores used with parallel
#' @return data.table object
#' @export
#' @docType methods
#' @rdname segvis-joinProfiles
#' @examples
#' \dontrun{ joinProfiles(segvis_object,bw = 151,mc=8)}
setGeneric("joinProfiles",
  function(object,bw,mc)
  standardGeneric("joinProfiles")
)

#' @title findSummit method for segvis class
#' @description This method finds the summit of each region whenever is possible, otherwise returns NA. The regions of the segvis object don't need to have the same width
#' @param object segvis object
#' @param bw, the bandwidth used to smooth the profiles, must be and odd number less or equal than the maxBandwidth
#' @param mc, the number of cores used with parallel
#' @return Numeric vector
#' @export
#' @docType methods
#' @rdname segvis-findSummit 
#' @examples
#' \dontrun{ findSummit(profile_object,bw = 151,mc=8)}
setGeneric("findSummit",
  function(object,bw,mc)
  standardGeneric("findSummit")
)

#' @title Create a segvis_block object
#' @description Constructor for segvis_block class, it takes both name and region from \code{object}
#' @param object segvis object
#' @param bw Numeric value with the bandwidth used to smooth the profiles, must be and odd number less or equal than the maxBandwidth
#' @param mc Numeric value with the number of cores used with parallel
#' @return a sevis_block object
#' @export
#' @docType methods
#' @rdname Segvis_block
#' @seealso \code{\link{segvis_block-class}}
#' @examples
#' \dontrun{Segvis_block(profile_object,bw = 151,mc=8)}
setGeneric("Segvis_block",
  function(object,bw,mc)
  standardGeneric("Segvis_block")
)           

#' @title Get and set method of reads objects
#' @description \code{readsF} operates on the forward reads, while \code{readsR} operates on the backward reads
#' @param object Either a \code{segvis} o \code{reads} object
#' @param value A data.table converted from a \code{GRanges} object
#' @export
#' @docType methods
#' @rdname methods-reads-gs
setGeneric("readsF",
  function(object)
  standardGeneric("readsF")
)           

#' @export
#' @docType methods
#' @rdname methods-reads-gs
setGeneric("readsF<-",
  function(object,value)
  standardGeneric("readsF<-")
)           

#' @export
#' @docType methods
#' @rdname methods-reads-gs
setGeneric("readsR",
  function(object)
  standardGeneric("readsR")           
)

#' @export
#' @docType methods
#' @rdname methods-reads-gs
setGeneric("readsR<-",
  function(object,value)
  standardGeneric("readsR<-")
)           

# Generic methods for segvis_block class

## Get methods

#' @title cover_table method for segvis_block class
#'
#' @description This methods determine how to get and set attributes for the segvis_block class. The structure is similar to the one use by bioconductor objects
#'
#' @param object An object of class segvis_block
#'
#' @return Returns the value with the same name as the function
#'
#' @export
#'
#' @docType methods
#'
#' @seealso \code{\link{segvis-class}}
#'
#' @rdname methods-segvis_block-gs
setGeneric("cover_table",
  function(object)
  standardGeneric("cover_table")           
)        

#' @param value The class of this parameter depends on the function, \code{name<-} is a string; for \code{bandwidth<-} and  \code{normConst} is a numeric value; for \code{regions<-} is a GRanges object and for \code{cover_table<-} is a data.table
#'
#' @return Returns a segvis object with the value indicated by the functions being updated with value
#'
#' @rdname methods-segvis_block-gs
setGeneric("cover_table<-",
  function(object,value)
  standardGeneric("cover_table<-")
)           

#' @export
#' @docType methods
#' @rdname methods-segvis_block-gs
setGeneric("bandwidth",
  function(object)
  standardGeneric("bandwidth")           
)           

#' @export
#' @docType methods
#' @rdname methods-segvis_block-gs
setGeneric("bandwidth<-",
  function(object,value)
  standardGeneric("bandwidth<-")
)           

#' @export
#' @docType methods
#' @rdname methods-segvis_block-gs
setGeneric("normConst",
  function(object)
  standardGeneric("normConst")           
)           

#' @export
#' @docType methods
#' @rdname methods-segvis_block-gs
setGeneric("normConst<-",
  function(object,value)
  standardGeneric("normConst<-")
)           

#' @title countReads method for segvis class
#' @description Counts the number of reads considered  in object
#' @param object segvis object
#' @return The number of reads in the bam file considered for the \code{segvis} object
#' @export
#' @docType methods
#' @seealso \code{\link{segvis}}
#' @rdname segvis-countReads
setGeneric("countReads",
  function(object)
  standardGeneric("countReads")
)           

#' @title summarize method for segvis_block class
#' @description This methods summarizes the coverage of the data when all the regions have the same width. To do so, it pilles up all tagCounts and applies FUN to each coordinate.
#' @param object segvis_block object
#' @param FUN Function used to summarize the tagCounts, for example if \code{FUN = mean}, then it will returns a vector with the mean of all tagCounts by centered coordinate. In particular, \code{FUN} needs to take vector argument + optional parameters
#' @param ... Rest of the arguments needed to use \code{FUN}
#' @return A numeric vector
#' @export
#' @docType methods
#' @seealso \code{\link{segvis_block-class}}
#' @rdname methods-segvis_block-summarize
setGeneric("summarize",
  function(object,FUN,...)
  standardGeneric("summarize"))


#' @title normalize method for segvis_block class
#' @description This method normalizes the segvis_block tagCounts by multiplying it times value.
#' It normalizes the tagCounts by multiplying it times base / normConst
#' @param object segvis_block object
#' @param value A positive number or a missing value. If is missing is going to normalize as if they were 1 million reads in the experiment
#' @param base A positive integer considering the base used to normalize the tagCounts
#' @return A segvis_block object where the cover_table contains normalized tagCounts
#' @export
#' @docType methods
#' @seealso \code{\link{segvis_block-class}}
#' @rdname methods-segvis_block-normalize
setGeneric("normalize",
  function(object,...)
  standardGeneric("normalize")
)           

#' @title subset method for segvis_block class
#' @description This method works similarly to the subset of IRanges, GenomicRanges, GenomicAlignments, etc. Althought it doesn't consider the select parameter.
#' @param object segvis_block object
#' @param condition This is an expression considering the characteristics taht the subset need to satisfy
#' @return Returns a segvis_block object with the same parameters as object except regions and profileMat which are filtered to satisfy the conditions on condition.
#' @export
#' @docType methods
#' @seealso \code{\link{segvis_block-class}}
#' @rdname subset
setGeneric("subset",
  function(object,condition)
  standardGeneric("subset"))

#' @title addColumn method for segvis_block class
#' @description This method helps to add a new column to the profile matrix object, is works similarly than using the $ operator over \code{regions(object)}
#' @param object segvis_block object
#' @param name Character with the name of the column to be add
#' @param col Vector of the same lengh as \code{regions(object)}
#' @return Returns a segvis_block object with the new column added to \code{regions(object)}
#' @export
#' @docType methods
#' @seealso \code{\link{segvis_block-class}}
#' @rdname methods-segvis_block-addColumn
setGeneric("addColumn",
  function(object,name,col)
  standardGeneric("addColumn"))           

#' @title plot.profile method for profileMatrixList class
#' @description This method returns the ggplot object used later to print the plot
#' @param object profileMatrixList object
#' @param condition The condition used to filter the profileMatrix objects, if there is no condition then it doesn't filter the objects
#' @param coord Numeric vector representing the coordinates used in the x-axis of the plot, by default considers the natural index \code{1:ncol(profileMatrix)}
#' @param trim Numeric. By default is set to zero. If trim =0, it calculated the row-wise mean, if trim >= 0.5 then is calculates the row-wise median, otherwise it calculate a row-wise trimmed mean.
#' @return The ggplot object used to print the plot
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrixList-class}}
#' @rdname plot.profiles
setGeneric("plot.profiles",
  function(object,...)
  standardGeneric("plot.profiles"))          

#' @title mergeList method for profileMatrixList class
#' @description This methods returns a profileMatrix class object which contains the regions and coverage curves of the list merged
#' @details In this case, we are going to merge by consider the priority implied by the order in list. Which means if there is an overlap between a region of object[[1]] with a region of object[[2]], the algorithm is going to choose the region of object[[1]] and so on.
#' @param object profileMatrixList with the profile matrices to be merged
#' @return A profileMatrix object with the regions and matrices merged
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrixList-class}}
#' @rdname mergeList
setGeneric("mergeList",
  function(object)
  standardGeneric("mergeList")
)           
