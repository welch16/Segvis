## Generic methods for SegvizData class

##' files method
##'
##' files returns a character vectors with the files used in the SegvizData
##' object.
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{files} method returns a character vector with the files
##' used in the \code{SegvizData} object.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname files-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' files(segviz)
##'
setGeneric("files",
      function(object)
      standardGeneric("files")
)

##' is_pet method
##'
##' is_pet returns a logical vector indicating if the aligned reads are paired
##' end
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{is_pet} method returns a logical vector indicating which
##' aligned reads files are paired ended.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname is_pet-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' is_pet(segviz)
##'
setGeneric("is_pet",
      function(object)
      standardGeneric("is_pet")
)

##' frag_len method
##'
##' frag_len returns a numeric vector with the fragment length used to extend
##' the 5' ends of the aligned reads.
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{frag_len} method returns a numeric vector with the fragment
##' length used to extend the 5' ends of the the aligned reads. If the reads
##' were paired ended, then a value of \code{NA} is assigned, but the parameter is
##' not used.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname frag_len-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' frag_len(segviz)
##'
setGeneric("frag_len",
        function(object)
        standardGeneric("frag_len")
)

##' covers method
##'
##' covers returns a list with the \code{coverage} of each file in the
##' \code{SegvizData} object.
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{covers} method returns a list of aggregated coverage for
##' all the aligned reads files considered.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname covers-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' covers(segviz)
##'
setGeneric("covers",
      function(object)
      standardGeneric("covers")
)

##' fwd_covers method
##'
##' fwd_covers returns a list with the \code{coverage} of each file in the
##' \code{SegvizData} object, calculated exclusively with the forward strand
##' reads.
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{fwd_covers} method returns a list of aggregated coverage
##' for all the aligned reads files considered, by using the reads in the forward
##' strand exclusively.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname fwd_covers-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' fwd_covers(segviz)
##'
setGeneric("fwd_covers",
      function(object)
      standardGeneric("fwd_covers")
)

##' bwd_covers method
##'
##' bwd_covers returns a list with the \code{coverage} of each file in the
##' \code{SegvizData} object, calculated exclusively with the backward strand
##' reads.
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{bwd_covers} method returns a list of aggregated coverage
##' for all the aligned reads files considered, by using the reads in the reverse
##' strand exclusively.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname bwd_covers-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' bwd_covers(segviz)
##'
setGeneric("bwd_covers",
        function(object)
        standardGeneric("bwd_covers")
)

##' nreads method
##'
##' nreads returns a numeric vector with the number of the depth of all the
##' experiments considered.
##'
##' @param object a \code{SegvizData} object.
##' @return The \code{nreads} method returns a numeric vector with the depth
##' of each experiment.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##'
##' @rdname nreads-methods
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' nreads(segviz)
##'
setGeneric("nreads",
      function(object)
      standardGeneric("nreads")
)

##' find_summits method
##'
##' find_summits returns a numeric vector with the same length of the
##' \code{GRanges} with the summits according to one of the files.
##'
##' @param object a \code{SegvizData} object.
##' @param which.file a numeric value indicating which file is going to be
##' used to calculate the summits. By default it takes the first file.
##' @param mc.cores a numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##'
##' @return The \code{find_summits} method returns a numeric vector with the
##' genomic coordinate with the summit in the \code{GRanges} according to one
##' of the files in object.
##'
##' @export
##' @rdname find_summits-methods
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' segviz$summit = find_summits(segviz,which.file = 1 ,mc.cores = 2)
##'
setGeneric("find_summits",
      function(object,...)
      standardGeneric("find_summits"))

##' overlap_matrix method
##'
##' overlap_matrix returns a matrix of zeros and ones indicating if the regions
##' in the \code{GRanges} object overlap a region in the respective bedfile
##'
##' @param object a \code{SegvizData} object.
##' @param bedfiles a character vector with the files in bed fortmat.
##' @param colnames a character vector with the names of the bed files. The
##' default value for this is \code{basename(bedfiles)}.
##'
##' @return The \code{overlap_matrix} method returns an overlap matrix where each
##' entry indicates if the row-region overlaps any of the intervals in the
##' bedfiles.
##'
##' @export
##' @rdname overlap_matrix-methods
##' @docType methods
##'
##' @seealso \code{\link{SegvizData-class}}
##' @examples
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' bedfiles = list.files(system.file("extdata","example",package = "Segvis"),
##'   pattern = "narrow",full.names = TRUE)[-1]
##'   overlap_matrix(segviz,bedfiles,c("dnase1","dnase2"))
##'
setGeneric("overlap_matrix",
      function(object,bedfiles,...)
      standardGeneric("overlap_matrix"))

##' DT_region method
##'
##' DT_region returns a \code{data.table} with the info. necessary to plot
##' the coverage of the aligned read files across the
##' \code{peak_id}th element of the \code{SegvizData} object.
##'
##' @param object a \code{SegvizData} object.
##' @param peak_id a numeric value indicating which region to plot.
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param normalize a logical value indicating if the profiles are normalized
##' respect to its number of reads.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##'
##' @return The \code{DT_region} method returns a \code{data.table} object with
##' the info. necessary to be plotted.
##' @rdname DT_region-methods
##'
##'@seealso \code{\link{SegvizData-class}},\code{\link{plot_region-methods}}
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' DT_region(segviz, peak_id = 1)
setGeneric("DT_region",
           function(object,peak_id,...)
             standardGeneric("DT_region"))

##' DT_profile method
##'
##' DT_profile return a \code{data.table} object with the summarized functional
##' profiles (summarized by using the argument \code{FUN}). For each coverage,
##' file and distance from the anchor coordinate it applies FUN to the tags
##' vector for all regions.
##'
##' @param object a \code{SegvizData} object.
##' @param FUN a function to summarize the profiles by coordinate. The default
##' value is the \code{mean}.
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##' \code{DT_profile} always normalizes the signal.
##' @param mc.cores a numeric value indicating the number of multi-cores to be
##' used.
##'
##' @return The \code{DT_profile} method returns a \code{data.table} object with
##' the info. necessary to be plotted.
##' @rdname DT_profile-methods
##'
##'@seealso \code{\link{SegvizData-class}},\code{\link{plot_profile-methods}}
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' DT_profile(segviz)
setGeneric("DT_profile",
           function(object,...)
             standardGeneric("DT_profile"))

##' plot_profile method
##'
##' plot_profile return a \code{ggplot2} object with the plot of the summarized
##' functional profiles (summarized by using the argument \code{FUN}). For each
##' coverage, file and distance from the anchor coordinate it applies FUN to
##' the tags vector for all regions.
##'
##' @param object a \code{SegvizData} object.
##' @param FUN a function to summarize the profiles by coordinate. The default
##' value is the \code{mean}.
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##' \code{DT_profile} always normalizes the signal.
##' @param mc.cores a numeric value indicating the number of multi-cores to be
##' used.
##'
##' @return The \code{plot_profile} method returns a \code{ggplot2} object with
##' the summarized functional profiles.
##' @rdname plot_profile-methods
##'
##'@seealso \code{\link{SegvizData-class}},\code{\link{DT_profile-methods}}
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' plot_profile(segviz)
setGeneric("plot_profile",
           function(object,...)
             standardGeneric("plot_profile"))

##' plot_heatmap method
##'
##' plot_heatmap plots a heatmap that summarizes the signal profiles of all files
##' in the \code{SegvizData} object as a function of the distance to the
##' center of the regions.
##'
##' @param object a \code{SegvizData} object.
##' @param which.cluster an integer indicating which signal is going to be used to
##' cluster the profiles. By default uses the first file loaded.
##' @param dist_method the distance measure to be used. This must be one of
##' "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
##' @param clust_method the agglomeration method to be used. This should be
##' one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
##' "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##' \code{DT_profile} always normalizes the signal.
##' @param mc.cores a numeric value indicating the number of multi-cores to be
##' used.
##'
##' @return The \code{plot_heatmap} plots a heatmap that summarizes the signal
##' profiles of all files in the \code{SegvizData} object.
##' @rdname plot_heatmap-methods
##'
##'@seealso \code{\link{SegvizData-class}}
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' plot_heatmap(segviz)
setGeneric("plot_heatmap",
           function(object,...)
             standardGeneric("plot_heatmap"))



##' plot_region method
##'
##' plot_region plots the coverage of the aligned read files across the
##' \code{peak_id}th element of the \code{SegvizData} object.
##'
##' @param object a \code{SegvizData} object.
##' @param peak_id a numeric value indicating which region to plot.
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param normalize a logical value indicating if the profiles are normalized
##' respect to its number of reads.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##'
##' @return The \code{plot_region} method returns a \code{ggplot} object with
##' the desired plot.
##' @rdname plot_region-methods
##'
##'@seealso \code{\link{SegvizData-class}}
##'
##' ## load SegvizData object
##' load(system.file("extdata","example","segviz.RData",package = "Segvis"))
##' plot_region(segviz, peak_id = 1)
setGeneric("plot_region",
        function(object,peak_id,...)
        standardGeneric("plot_region"))

############## possibly methods to be removed

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

