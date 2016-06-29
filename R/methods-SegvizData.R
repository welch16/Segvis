

##' @rdname files-methods
##' @aliases files
##' @docType methods
##' @exportMethod files
setMethod("files",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@files)

##' @rdname is_pet-methods
##' @aliases is_pet
##' @docType methods
##' @exportMethod is_pet
setMethod("is_pet",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@is_pet)

##' @rdname frag_len-methods
##' @aliases frag_len
##' @docType methods
##' @exportMethod frag_len
setMethod("frag_len",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@frag_len)

##' @rdname covers-methods
##' @aliases covers
##' @docType methods
##' @exportMethod covers
setMethod("covers",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@covers)

##' @rdname fwd_covers-methods
##' @aliases fwd_covers
##' @docType methods
##' @exportMethod fwd_covers
setMethod("fwd_covers",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@fwd_covers)

##' @rdname bwd_covers-methods
##' @aliases bwd_covers
##' @docType methods
##' @exportMethod bwd_covers
setMethod("bwd_covers",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@bwd_covers)

##' @rdname nreads-methods
##' @aliases nreads
##' @docType methods
##' @exportMethod nreads
setMethod("nreads",
      signature = signature(object = "SegvizData"),
      definition = function(object)object@nreads)

##' @rdname find_summits-methods
##' @aliases find_summits
##' @docType methods
##' @exportMethod find_summits
setMethod("find_summits",
      signature = signature(object = "SegvizData"),
      definition = function(object,which.file = 1,
          mc.cores = getOption("mc.cores",2L)){
        stopifnot(is.numeric(which.file),which.file >= 1,
                  which.file <= length(files(object)))
        cover = covers(object)[[which.file]]
        subs = cover[object]
        sm = mclapply(subs, .rle_summit,mc.cores = mc.cores)
        names(sm) = NULL
        start(object) + do.call(c,sm)
      })

##' @rdname overlap_matrix-methods
##' @aliases overlap_matrix
##' @docType methods
##' @exportMethod overlap_matrix
setMethod("overlap_matrix",
      signature = signature(object = "SegvizData",bedfiles = "character"),
      definition = function(object,bedfiles,
                            colnames = basename(bedfiles)){
        stopifnot(is.character(bedfiles),all(file.exists(bedfiles)))

        regions = lapply(bedfiles,readBedFile)
        overlaps = lapply(regions,function(x)overlapsAny(object,x))

        mat = DataFrame(ifelse(do.call(cbind,overlaps),1,0))
        colnames(mat) = colnames

        mat

      })




