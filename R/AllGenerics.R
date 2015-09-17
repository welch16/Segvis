
# generic methods for various classes

##' name methods
##'
##' name returns a string with the name of the segvis object
##' 
##' @param object A \code{segvis} object
##'
##' @return A string with the name of the object
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}
##' @rdname name-methods
##' @examples
##' \dontrun{
##'
##' name(segvis)
##' name(segvis) <- "my_regions"
##' 
##' }
setGeneric("name",
  function(object)
  standardGeneric("name")
)

##' name<- assigns a new a name to the segvis object
##'
##' @param value A character with the name of the object
##'
##' @return A segvis object
##' @rdname name-methods
setGeneric("name<-",
  function(object,value)
  standardGeneric("name<-")
)           

##' regions methods
##'
##' regions return a GenomicRanges object with the regions used by segvis
##' 
##' @param object A \code{segvis} object
##'
##' @return A GenomicRanges object
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}
##' @rdname regions-methods
##' @examples
##' \dontrun{
##'
##' regions(segvis)
##' regions(segvis) <- gr  ## gr is a GenomicRanges object
##' 
##' }
setGeneric("regions",
  function(object)
  standardGeneric("regions")
)

##' regions<- assigns a GenomicRanges object as the regions used by segvis
##'
##' @param value A GenomicRanges object
##'
##' @return A segvis object
##' @rdname regions-methods
setGeneric("regions<-",
  function(object,value)
  standardGeneric("regions<-")
)

##' file methods
##'
##' file return the file name used in the segvis object
##' 
##' @param object A \code{segvis} object
##'
##' @return A string
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}
##' @rdname file-methods
##' @examples
##' \dontrun{
##'
##' file(segvis)
##' file(segvis) <- "path/to/file/myFile.bam"
##' 
##' }
setGeneric("file",           
  function(object)
  standardGeneric("file")           
)


##' file<- assigns a string with the file name of the reads to be used 
##'
##' @param value A character value with the name of bam the file with the reads 
##'
##' @return A segvis object
##' @rdname file-methods
setGeneric("file<-",
  function(object,value)
  standardGeneric("file<-")
)

##' maxBandwidth methods
##'
##' maxBandwidth returns a numeric value representing the possible bandwidths upper bound to smooth the coverage profiles.
##' 
##' @param object A \code{segvis} object
##'
##' @return A integer value
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}
##' @rdname maxBandwidth-methods
##' @examples
##' \dontrun{
##'
##' maxBandwidth(segvis)
##' maxBandwidth(segvis) <- 201
##' 
##' }
setGeneric("maxBandwidth",
  function(object)
  standardGeneric("maxBandwidth")           
)


##' maxBandwidth<- assigns an integer value the bandwidth upper bound to use when smoothing the profiles in segvis_block object
##'
##' @param value A numeric value representing the maximum bandwidth, it need to be an odd number
##'
##' @return A segvis object
##' @rdname maxBandwidth-methods
setGeneric("maxBandwidth<-",
  function(object,value)
  standardGeneric("maxBandwidth<-")
)           

##' fragLen methods
##'
##' fragLen returns a numeric value representing the fragment length used to extend the fragment reads in the Single Ended case (i.e. isPET = FALSE)
##' 
##' @param object A \code{segvis} object
##'
##' @return A integer value
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}
##' @rdname fragLen-methods
##' @examples
##' \dontrun{
##'
##' fragLen(segvis)
##' fragLen(segvis) <- 200
##' 
##' }
setGeneric("fragLen",
  function(object)
  standardGeneric("fragLen")           
)

##' fragLen<- assigns an integer value representing the bp to extend the fragments in SET case
##'
##' @param value A numeric value representing the fragment length
##'
##' @return A segvis object
##' @rdname fragLen-methods
setGeneric("fragLen<-",
  function(object,value)
  standardGeneric("fragLen<-")
)           

##' chr methods
##'
##' chr returns a character vector naming the chromosomes used in the analysis
##'
##' @param object A \code{segvis} object
##'
##' @return A character vector
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}, \code{\link{buildSegvis}}
##' @rdname chr-methods
##' @examples
##' \dontrun{
##'
##' chr(segvis)
##' chr(segvis) <- c("chr1","chr2",...,"chrX","chrY") ## this is not correct actually
##' 
##' }
setGeneric("chr",
  function(object)
  standardGeneric("chr")
)


##' chr<- assigns a logical flag to the segvis object
##'
##' @param value A character vector with the chromosomes to use in the analysis
##'
##' @return A segvis object
##' @rdname chr-methods
setGeneric("chr<-",
  function(object,value)
  standardGeneric("chr<-")           
)


##' isPET methods
##'
##' isPET returns a logical flag wheter the reads to consider are Paired Ended or not
##'
##' @param object A \code{segvis} object
##'
##' @return A logical flag
##'
##' @export
##' @docType methods
##' @seealso \code{\link{segvis-class}}
##' @rdname isPET-methods
##' @examples
##' \dontrun{
##'
##' isPET(segvis)
##' isPET(segvis) <- TRUE
##' 
##' }
setGeneric("isPET",
  function(object)
  standardGeneric("isPET")           
)
           
##' isPET<- assigns a logical flag to the segvis object
##'
##' @param value A logical flag that represents wheter the reads to use are Paired Ended or not
##'
##' @return A segvis object
##' @rdname isPET-methods
setGeneric("isPET<-",
  function(object,value)
  standardGeneric("isPET<-")           
)           




##' profiles methods for segvis class
##'
##' Returns the coverage for each region considered
##'
##' @param object A segvis object
##'
##' @return RleList. A list made of an Rle object for each region
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis-class}} and \code{\link{getCoverage}}
##'
##' @rdname profiles-methods
##'
##' @examples
##' \dontrun{
##'
##' profiles(segvis)
##'
##' }
setGeneric("profiles",
  function(object)
  standardGeneric("profiles")
)  
    
##' loadReads method for segvis class
##'
##' Load the fragment stored in the file slot of the segvis object. The reads are divided by chromosome and by strand.
##' @param object segvis object
##'
##' @param mc numeric, the number of cores used with parallel
##'
##' @return segvis object
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{reads-class}}
##'
##' @rdname loadReads-methods
##'
##' @examples
##' \dontrun{
##'
##' segvis <- loadReads(segvis,mc=8)
##'
##' }
setGeneric("loadReads",
  function(object,mc)
  standardGeneric("loadReads")
)           

##' matchReads methods for segvis class
##'
##' Match the reads to the extended regions stored in the object slots.
##'
##' @param object segvis object
##'
##' @param mc numeric, the number of cores used with parallel
##'
##' @return segvis object
##' 
##' @export
##' @docType methods
##'
##' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{reads-class}}
##'
##' @rdname matchReads-methods
##'
##' @examples
##' \dontrun{
##'
##' matchReads(segvis_object,mc=8)
##'
##' }
setGeneric("matchReads",
  function(object,mc)
  standardGeneric("matchReads")
)           


##' getCoverage method for profile class
##'
##' Calculate  the coverage using the reads matched by the matchReads method. It returns a Rle object for each region.
##'
##' @param object segvis object
##' 
##' @param mc numeric, the number of cores used with parallel
##'
##' @return segvis object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{loadReads}}, \code{\link{matchReads}} and \code{\link{findSummit}}
##' 
##' @rdname getCoverage-methods
##'
##' @examples
##' \dontrun{
##' 
##' getCoverage(segvis_object,mc=8)
##'
##' }
setGeneric("getCoverage",
  function(object,mc)
  standardGeneric("getCoverage")           
)


##' joinProfiles method for segvis class
##'
##' Joins all the coverages into a data.table with columns chr|match|coord|tagCounts
##'
##' @param object segvis object
##'
##' @param bw Numeric value with the bandwidth used to smooth the coverage profiles, must be and odd number less or equal than the maxBandwidth
##'
##' @param mc Numeric value with the number of cores used with parallel
##'
##' @return A data.table object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}},\code{\link{segvis-class}}
##'
##' @rdname joinProfiles-methods
##'
##' @examples
##' \dontrun{
##'
##' joinProfiles(segvis_object,bw = 151,mc=4)
##' 
##' }
setGeneric("joinProfiles",
  function(object,bw,mc)
  standardGeneric("joinProfiles")
)

##' findSummit method for segvis class
##' 
##' This method finds the summit of each region whenever is possible, otherwise returns NA. The regions of the segvis object don't need to have the same width
##' 
##' @param object segvis object
##' 
##' @param bw, the bandwidth used to smooth the profiles, must be and odd number less or equal than the maxBandwidth
##' 
##' @param mc, the number of cores used with parallel
##'
##' @return Numeric vector with the estimated summits
##' 
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{getCoverage}},\code{\link{segvis-class}}
##' @rdname findSummit-methods
##'
##' @examples
##' \dontrun{
##' 
##' summits <- findSummit(profile_object,bw = 151,mc=8)
##' 
##' }
setGeneric("findSummit",
  function(object,bw,mc)
  standardGeneric("findSummit")
)

##' Segvis_block constructor from a segvis object
##'
##' This method returns a segvis_block object built out of the arguments in a segvis object
##'
##' @param object segvis object
##'
##' @param bw Numeric value with the bandwidth used to smooth the coverage profiles, must be and odd number less or equal than the maxBandwidth
##'
##' @param mc Numeric value with the number of cores used with parallel
##'
##' @return a segvis_block object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}},\code{\link{segvis-class}}
##'
##' @rdname Segvis_block-methods
##'
##' @examples
##' \dontrun{
##'
##' mc <- 4
##' segvis_block <- Segvis_block(segvis,101,mc)
##'
##' }
setGeneric("Segvis_block",
  function(object,bw,mc)
  standardGeneric("Segvis_block")
)           

##' readsF methods
##'
##' readsF returns a list of forward reads of segvis or reads objects
##'
##' @param object Either a \code{segvis} or a \code{reads} object
##'
##' @return A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @export
##' @docType methods
##' @seealso \code{\link{reads-class}},\code{\link{segvis-class}}
##' @rdname readsF-methods
##' @examples
##' \dontrun{
##' readsF(segvis)
##' readsF(reads)
##' readsF(segvis) <- new_reads
##' }
setGeneric("readsF",
  function(object)
  standardGeneric("readsF")
)           

##' readsF<- assisgn a list of forward reads to a reads or to a segvis object
##'
##' @param value A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @return A reads or a segvis object with modified forward reads
##' @rdname readsF-methods
setGeneric("readsF<-",
  function(object,value)
  standardGeneric("readsF<-")
)           


##' readsR methods
##'
##' readsR returns a list of backward reads of segvis or reads objects
##'
##' @param object Either a \code{segvis} o \code{reads} object
##'
##' @return A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @export
##' @docType methods
##' @seealso \code{\link{reads-class}},\code{\link{segvis-class}}
##' @rdname readsR-methods
##' @examples
##' \dontrun{
##' readsR(segvis)
##' readsR(reads)
##' readsR(segvis) <- new_reads
##' }
setGeneric("readsR",
  function(object)
  standardGeneric("readsR")           
)

##' readsR<- assisgn a list of backward reads to a reads or to a segvis object
##'
##' @param value A list of data.table with same format as a GenomicRanges object with an additional match column
##'
##' @return A reads or a segvis object with modified backward reads
##' @rdname readsR-methods
setGeneric("readsR<-",
  function(object,value)
  standardGeneric("readsR<-")
)           

# Generic methods for segvis_block class
# Get methods

##' cover_table methods for segvis_block class
##'
##' cover_table method returns coverage table of a segvis_block object
##'
##' @param object segvis_block object
##'
##' @return The coverage data.table of a segvis_block object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname cover_table-methods
##'
setGeneric("cover_table",
  function(object)
  standardGeneric("cover_table")           
)        

##' cover_table<- assigns a new coverage table to a segvis_block object
##'
##' @param value data.table with the new coverage of the se
##'
##' @return A segvis_block object with the bandwidth replcated
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname cover_table-methods
##'
setGeneric("cover_table<-",
  function(object,value)
  standardGeneric("cover_table<-")
)


##' bandwidth methods for segvis_block class
##'
##' bandwidth method returns the smoothing bandwidth of a segvis_block object
##'
##' @param object segvis_block object
##'
##' @return The smoothing bandwidth of the segvis_block object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname bandwidth-methods
##'
##' @examples
##' \dontrun{
##'
##' bandwidth(segvis_block)
##' bandwidth(segvis_block) <- 101
##'
##' }
setGeneric("bandwidth",
  function(object)
  standardGeneric("bandwidth")           
)           


##' bandwidth<- assigns a new smoothing bandwidth to a segvis_block object
##'
##' @param value Numeric value
##'
##' @return A segvis_block object with the bandwidth replcated
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname bandwidth-methods
##'
setGeneric("bandwidth<-",
  function(object,value)
  standardGeneric("bandwidth<-")
)           

##' normConst methods for segvis_block class
##'
##' normConst method returns the normalizing constant of a segvis_block object
##'
##' @param object segvis_block object
##'
##' @return The normalizing constant of the segvis_block object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname normConst-methods
##'
##' @examples
##' \dontrun{
##'
##' normConst(segvis_block)
##' normConst(segvis_block) <- 100
##' normConst(segvis_block) <- countReads(segvis)
##'
##' }
setGeneric("normConst",
  function(object)
  standardGeneric("normConst")           
)           


##' normConst<- assigns a new normalizing constant to a segvis_block object
##'
##' @param value Numeric value
##'
##' @return A segvis_block object with the normConst replcated
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname normConst-methods
##'
setGeneric("normConst<-",
  function(object,value)
  standardGeneric("normConst<-")
)           

##' countReads method for segvis class
##'
##' Counts the number of reads considered  in object
##'
##' @param object segvis object
##'
##' @return The number of reads in the bam file considered for the \code{segvis} object
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis-class}}
##'
##' @rdname countReads-methods
##'
##' @examples
##' \dontrun{
##'
##' countReads(segvis)
##'
##' }
setGeneric("countReads",
  function(object)
  standardGeneric("countReads")
)           

##' summarize method for segvis_block class
##'
##' This methods summarizes the coverage of the data when all the regions have the same width. To do so, it pilles up all tagCounts and applies FUN to each coordinate.
##'
##' @param object segvis_block object
##'
##' @param FUN Function used to summarize the tagCounts, for example if \code{FUN = mean}, then it will returns a vector with the mean of all tagCounts by centered coordinate. In particular, \code{FUN} needs to take vector argument + optional parameters
##'
##' @param ... Rest of the arguments needed to use \code{FUN}
##'
##' @return A numeric vector
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname summarize-methods
##'
##' @examples
##' \dontrun{
##'
##' summarize(segvis_block,mean)
##' summarize(segvis_block,mean,trim=.25)
##' summarize(segvis_block,median)
##'
##' }
setGeneric("summarize",
  function(object,FUN,...)
  standardGeneric("summarize"))


##' normalize method for segvis_block class
##'
##' This method normalizes the segvis_block tagCounts by multiplying it times value. It normalizes the tagCounts by multiplying it times base / normConst
##'
##' @param object segvis_block object
##'
##' @param base Numeric value, indicating to which scale the profiles are gonna be scaled
##' 
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname normalize-methods
##'
##' @examples
##' \dontrun{
##'
##' normalize(segvis_block)
##' normalize(segvis_block,base = 1e9) ## rpkm
##'
##' }
setGeneric("normalize",
  function(object,base)
  standardGeneric("normalize")
)           

##' subset_block method for segvis_block class
##'
##' This method works similarly to the subset of IRanges, GenomicRanges, GenomicAlignments, etc. Althought it doesn't consider the select parameter.
##'
##' @param object segvis_block object
##'
##' @param condition This is an expression considering the characteristics taht the subset need to satisfy
##'
##' @return Returns a segvis_block object with the same parameters as object except regions and cover_table which are filtered to satisfy the conditions on condition.
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname subset_block-methods
##'
##' @examples
##' \dontrun{
##'
##' subset_block(segvis_block, cond == TRUE)
##'
##' }
setGeneric("subset_block",
  function(object,condition)
  standardGeneric("subset_block"))

##' addColumn method for segvis_block class
##'
##' This method helps to add a new column to the profile matrix object, is works similarly than using the $ operator over \code{regions(object)}
##'
##' @param object segvis_block object
##'
##' @param name Character with the name of the column to be add
##'
##' @param col Vector of the same lengh as \code{regions(object)}
##'
##' @return Returns a segvis_block object with the new column added to \code{regions(object)}
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block-class}}
##'
##' @rdname addColumn-methods
##'
##' @examples
##' \dontrun{
##'
##' addColumn(segvis_block,"cond_name",cond_vec)
##'
##' }
setGeneric("addColumn",
  function(object,name,col)
  standardGeneric("addColumn"))           

##' plot_profiles method for segvis_block_list class
##'
##' This method returns the ggplot object used later to print the plot
##'
##' @param object segvis_block_list object
##'
##' @param condition The condition used to filter the profileMatrix objects, if there is no condition then it doesn't filter the objects
##'
##' @param coord Numeric vector representing the coordinates used in the x-axis of the plot, by default considers the natural index \code{1:ncol(segvis_block)}
##'
##' @param FUN Function used to summarize the profiles
##'
##' @param ... Additional arguments of \code{FUN}
##'
##' @param mc Number of cores used by parallel
##'
##' @return The ggplot object used to print the plot
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block_list-class}},\code{\link{plot_data}}
##'
##' @rdname plot_profiles-methods
##'
##' @examples
##' \dontrun{
##'
##'  plot_profiles(segvis_block_list, mean, trim = .05,mc = 8)
##'
##'  ## if the segvis_block_list have a column called row
##'  plot_profiles(segvis_block_list, FUN = function(x)x,condition = row == 1,mc = 8)
##'
##' }
setGeneric(name="plot_profiles",
  def=function(object,condition,coord,FUN,...,mc)
  standardGeneric("plot_profiles"))          

##' plot_data method for segvis_block_list class
##'
##' This method returns the data.table used to later make a ggplot. The idea is to used this table to make more complicated plots
##'
##' @param object segvis_block_list object
##'
##' @param condition The condition used to filter the profileMatrix objects, if there is no condition then it doesn't filter the objects
##'
##' @param coord Numeric vector representing the coordinates used in the x-axis of the plot, by default considers the natural index \code{1:ncol(segvis_block)}
##'
##' @param FUN Function used to summarize the profiles
##'
##' @param ... Additional arguments of \code{FUN}
##'
##' @param mc Number of cores used by parallel
##'
##' @return A data.table with the data used to make a plot
##'
##' @export
##'
##' @docType methods
##'
##' @seealso \code{\link{segvis_block_list-class}},\code{\link{plot_profiles}}
##'
##' @rdname plot_data-methods
##'
##' @examples
##' \dontrun{
##'
##'  plot_data(segvis_block_list, mean, trim = .05,mc = 8)
##'
##'  ## if the segvis_block_list have a column called row
##'  plot_data(segvis_block_list, FUN = function(x)x,condition = row == 1,mc = 8)
##'
##' }
setGeneric(name="plot_data",
  def=function(object,condition,coord,FUN,...,mc)standardGeneric("plot_data"))
