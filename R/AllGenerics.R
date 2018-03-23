## Generic methods for SegvisData class

##' files method
##'
##' files returns a character vectors with the files used in the SegvisData
##' object.
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{files} method returns a character vector with the files
##' used in the \code{SegvisData} object.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname files-methods
##
setGeneric("files",
      function(object)
      standardGeneric("files")
)

##' is_pet method
##'
##' is_pet returns a logical vector indicating if the aligned reads are paired
##' end
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{is_pet} method returns a logical vector indicating which
##' aligned reads files are paired ended.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname is_pet-methods
##
setGeneric("is_pet",
      function(object)
      standardGeneric("is_pet")
)

##' frag_len method
##'
##' frag_len returns a numeric vector with the fragment length used to extend
##' the 5' ends of the aligned reads.
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{frag_len} method returns a numeric vector with the fragment
##' length used to extend the 5' ends of the the aligned reads. If the reads
##' were paired ended, then a value of \code{NA} is assigned, but the parameter is
##' not used.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname frag_len-methods
##
setGeneric("frag_len",
        function(object)
        standardGeneric("frag_len")
)

##' covers method
##'
##' covers returns a list with the \code{coverage} of each file in the
##' \code{SegvisData} object.
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{covers} method returns a list of aggregated coverage for
##' all the aligned reads files considered.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname covers-methods
##
setGeneric("covers",
      function(object)
      standardGeneric("covers")
)

##' fwd_covers method
##'
##' fwd_covers returns a list with the \code{coverage} of each file in the
##' \code{SegvisData} object, calculated exclusively with the forward strand
##' reads.
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{fwd_covers} method returns a list of aggregated coverage
##' for all the aligned reads files considered, by using the reads in the forward
##' strand exclusively.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname fwd_covers-methods
##
setGeneric("fwd_covers",
      function(object)
      standardGeneric("fwd_covers")
)

##' bwd_covers method
##'
##' bwd_covers returns a list with the \code{coverage} of each file in the
##' \code{SegvisData} object, calculated exclusively with the backward strand
##' reads.
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{bwd_covers} method returns a list of aggregated coverage
##' for all the aligned reads files considered, by using the reads in the reverse
##' strand exclusively.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname bwd_covers-methods
##
setGeneric("bwd_covers",
        function(object)
        standardGeneric("bwd_covers")
)

##' nreads method
##'
##' nreads returns a numeric vector with the number of the depth of all the
##' experiments considered.
##'
##' @param object a \code{SegvisData} object.
##' @return The \code{nreads} method returns a numeric vector with the depth
##' of each experiment.
##'
##' @export
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##'
##' @rdname nreads-methods
##
setGeneric("nreads",
      function(object)
      standardGeneric("nreads")
)

##' find_summits method
##'
##' find_summits returns a numeric vector with the same length of the
##' \code{GRanges} with the summits according to one of the files.
##'
##' @param object a \code{SegvisData} object.
##' @param which.file a numeric value indicating which file is going to be
##' used to calculate the summits. By default it takes the first file.
##' @param mc.cores a numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##' @param ... Other parameters that may be used by \code{find_summit}.
##'
##' @return The \code{find_summits} method returns a numeric vector with the
##' genomic coordinate with the summit in the \code{GRanges} according to one
##' of the files in object.
##'
##' @export
##' @rdname find_summits-methods
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##
setGeneric("find_summits",
      function(object,...)
      standardGeneric("find_summits"))

##' overlap_matrix method
##'
##' overlap_matrix returns a matrix of zeros and ones indicating if the regions
##' in the \code{GRanges} object overlap a region in the respective bedfile
##'
##' @param object a \code{SegvisData} object.
##' @param bedfiles a character vector with the files in bed fortmat.
##' @param colnames a character vector with the names of the bed files. The
##' default value for this is \code{basename(bedfiles)}.
##' @param ... Any other additional parameters that \code{overlap_matrix} may need.
##'
##' @return The \code{overlap_matrix} method returns an overlap matrix where each
##' entry indicates if the row-region overlaps any of the intervals in the
##' bedfiles.
##'
##' @export
##' @rdname overlap_matrix-methods
##' @docType methods
##'
##' @seealso \code{\link{SegvisData-class}}
##
setGeneric("overlap_matrix",
      function(object,bedfiles,...)
      standardGeneric("overlap_matrix"))

##' DT_region method
##'
##' DT_region returns a \code{data.table} with the info. necessary to plot
##' the coverage of the aligned read files across the
##' \code{peak_id}th element of the \code{SegvisData} object.
##'
##' @param object a \code{SegvisData} object.
##' @param region a numeric value indicating which region to plot or a
##' \code{GRanges} object indicating a whole region to plot.
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param normalize a logical value indicating if the profiles are normalized
##' respect to its number of reads.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##' @param ... Any other additional parameters that \code{DT_region} may need.
##'
##' @return The \code{DT_region} method returns a \code{data.table} object with
##' the info. necessary to be plotted.
##' @rdname DT_region-methods
##'
##'@seealso \code{\link{SegvisData-class}}
##'
setGeneric("DT_region",
           function(object,region,...)
             standardGeneric("DT_region"))

##' DT_profile method
##'
##' DT_profile return a \code{data.table} object with the summarized functional
##' profiles (summarized by using the argument \code{FUN}). For each coverage,
##' file and distance from the anchor coordinate it applies FUN to the tags
##' vector for all regions.
##'
##' @param object a \code{SegvisData} object.
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
##' @param ... Any other additional parameters that \code{DT_profile} may need.
##'
##' @return The \code{DT_profile} method returns a \code{data.table} object with
##' the info. necessary to be plotted.
##' @rdname DT_profile-methods
##'
##'@seealso \code{\link{SegvisData-class}}
##'
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
##' @param object a \code{SegvisData} object.
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
##' @param ... Any other additional parameters that \code{plot_profile} may need.
##'
##' @return The \code{plot_profile} method returns a \code{ggplot2} object with
##' the summarized functional profiles.
##' @rdname plot_profile-methods
##'
##'@seealso \code{\link{SegvisData-class}}
##'
setGeneric("plot_profile",
           function(object,...)
             standardGeneric("plot_profile"))

##' plot_heatmap method
##'
##' plot_heatmap plots a heatmap that summarizes the signal profiles of all files
##' in the \code{SegvisData} object as a function of the distance to the
##' center of the regions.
##'
##' @param object a \code{SegvisData} object.
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
##' @param ... Any other additional parameters that \code{plot_heatmap} may need.
##'
##' @return The \code{plot_heatmap} plots a heatmap that summarizes the signal
##' profiles of all files in the \code{SegvisData} object.
##' @rdname plot_heatmap-methods
##'
##' @seealso \code{\link{SegvisData-class}}
##' 
setGeneric("plot_heatmap",
           function(object,...)
             standardGeneric("plot_heatmap"))



##' plot_region method
##'
##' plot_region plots the coverage of the aligned read files across the
##' \code{peak_id}th element of the \code{SegvisData} object.
##'
##' @param object a \code{SegvisData} object.
##' @param region a numeric value indicating which region to plot in the
##' \code{SegvisData} object or a \code{GRanges} object with a genomic region
##' that is not recorded in the \code{SegvisData} object
##' @param nameFiles a character vector with the shortened names that are going
##' to be used in the plot.
##' @param type a character value indicating which strand to use. By default
##' uses the aggregate coverage between both strands.
##' @param normalize a logical value indicating if the profiles are normalized
##' respect to its number of reads.
##' @param base a numeric value indicating the number of aligned reads to which
##' the signal is going to be normalized. The default value is \code{1e6}.
##' @param ... Any other additional parameters that \code{plot_region} may need.
##' @return The \code{plot_region} method returns a \code{ggplot} object with
##' the desired plot.
##' @rdname plot_region-methods
##'
##'@seealso \code{\link{SegvisData-class}}
##'
setGeneric("plot_region",
        function(object,region,...)
        standardGeneric("plot_region"))

