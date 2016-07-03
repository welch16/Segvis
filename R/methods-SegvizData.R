

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

##' @rdname DT_region-methods
##' @aliases DT_region
##' @docType methods
##' @exportMethod DT_region
setMethod("DT_region",
          signature = signature(object = "SegvizData",peak_id = "numeric"),
          definition = function(object,peak_id,
                                nameFiles = basename(files(object)),
                                type = "aggr",normalize = TRUE,
                                base = 1e6){
            if(is.numeric(peak_id) & !is.integer(peak_id)){
              peak_id = as.integer(peak_id)
            }

            stopifnot(is.character(type),
                      tolower(type) %in% c("aggr","fwd","bwd","both","all"))
            stopifnot(is.integer(peak_id),length(peak_id) == 1,
                      is.logical(normalize))
            stopifnot(1 <= peak_id , peak_id <= length(object))

            if(is.numeric(nameFiles))nameFiles = as.character(nameFiles)

            reg = object[peak_id]

            aggr_dt = NULL
            fwd_dt = NULL
            bwd_dt = NULL

            if(tolower(type) %in% c("aggr","all")){
              aggr_dt = mapply(.dt_cover,covers(object),
                               nreads(object),nameFiles,
                               MoreArgs = list(region = reg,
                                               st = "aggr",
                                               normalize = normalize,
                                               base = base),
                               SIMPLIFY = FALSE)
            }

            if(tolower(type) %in% c("fwd","both","all")){
              fwd_dt = mapply(.dt_cover,fwd_covers(object),
                              nreads(object),nameFiles,
                              MoreArgs = list(region = reg,
                                              st = "fwd",
                                              normalize = normalize,
                                              base = base),
                              SIMPLIFY = FALSE)
            }
            if(tolower(type) %in% c("bwd","both","all")){
              bwd_dt = mapply(.dt_cover,bwd_covers(object),
                              nreads(object),nameFiles,
                              MoreArgs = list(region = reg,
                                              st = "bwd",
                                              normalize = normalize,
                                              base = base),
                              SIMPLIFY = FALSE)
            }

            DT = rbindlist(c(fwd_dt,bwd_dt,aggr_dt))
            DT
          })

##' @rdname DT_profile-methods
##' @aliases DT_profile
##' @docType methods
##' @exportMethod DT_profile
setMethod("DT_profile",
          signature = signature(object = "SegvizData"),
          definition = function(object,FUN = mean,
                                nameFiles = basename(files(object)),
                                type = "aggr",
                                base = 1e6,
                                mc.cores = getOption("mc.cores",2L),...){

            stopifnot(is.character(type),
                      tolower(type) %in% c("aggr","fwd","bwd","both","all"))
            stopifnot(length(unique(width(object))) == 1)

            if(is.numeric(nameFiles))nameFiles = as.character(nameFiles)

            aggr_dt = NULL
            fwd_dt = NULL
            bwd_dt = NULL

            width = unique(width(object))
            ll = (width - 1) / 2

            if(tolower(type) %in% c("aggr","all")){
              aggr_dt = mapply(.dt_profile,covers(object),
                               nreads(object),nameFiles,
                               MoreArgs = list(regions = object,
                                               st = "aggr",
                                               base = base,
                                               len = ll,
                                               mc.cores = mc.cores),
                               SIMPLIFY = FALSE)
            }

            if(tolower(type) %in% c("fwd","both","all")){
              fwd_dt = mapply(.dt_profile,fwd_covers(object),
                               nreads(object),nameFiles,
                               MoreArgs = list(regions = object,
                                               st = "fwd",
                                               base = base,
                                               len = ll,
                                               mc.cores = mc.cores),
                               SIMPLIFY = FALSE)

            }
            if(tolower(type) %in% c("bwd","both","all")){

              bwd_dt = mapply(.dt_profile,bwd_covers(object),
                               nreads(object),nameFiles,
                               MoreArgs = list(regions = object,
                                               st = "bwd",
                                               base = base,
                                               len = ll,
                                               mc.cores = mc.cores),
                               SIMPLIFY = FALSE)

            }

            DT = rbindlist(c(fwd_dt,bwd_dt,aggr_dt))
            DT = DT[,FUN(tags,...),by = .(coord,name,type)]
            setnames(DT,names(DT),c("coord","name","type","tags"))
            DT
           })

##' @rdname plot_profile-methods
##' @aliases plot_profile
##' @docType methods
##' @exportMethod plot_profile
setMethod("plot_profile",
          signature = signature(object = "SegvizData"),
          definition = function(object,FUN = mean,
                                nameFiles = basename(files(object)),
                                type = "aggr",
                                base = 1e6,
                                mc.cores = getOption("mc.cores",2L),
                                ...){

            DT = DT_profile(object,FUN,nameFiles,type,base,mc.cores)
            pal = c("black",brewer.pal(9,"Set1"))

            if(tolower(type ) %in% c("aggr","fwd","bwd")){
              out = ggplot(DT,aes_string(x="coord",y= "tags",
                                         linetype = "name"))+
                geom_line()+xlab("Genomic Coordinates")+
                ylab("Normalized Signal")
            }else if(tolower(type) == "both"){
              pal = pal[-1]
              out = ggplot(DT,aes_string(x="coord",y= "tags",
                                         linetype = "name",
                                         colour = "type"))+
                geom_line()+xlab("Genomic Coordinates")+
                ylab("Normalized Signal")+
                scale_color_manual(values = pal)
            }else{

              out = ggplot(DT,aes_string(x="coord",y= "tags",
                                         linetype = "name",
                                         colour = "type"))+
                geom_line()+xlab("Genomic Coordinates")+
                ylab("Normalized Signal")+
                scale_color_manual(values = pal)
            }
            out

          })

##' @rdname plot_region-methods
##' @aliases plot_region
##' @docType methods
##' @exportMethod plot_region
setMethod("plot_region",
      signature = signature(object = "SegvizData",peak_id = "numeric"),
      definition = function(object,peak_id,
                            nameFiles = basename(files(object)),
                            type = "aggr",normalize = TRUE,
                            base = 1e6){
        DT = DT_region(object,peak_id,nameFiles = nameFiles,
                      type = type,normalize = normalize,
                      base = base)
        pal = c("black",brewer.pal(9,"Set1"))

        if(tolower(type ) %in% c("aggr","fwd","bwd")){
          out = ggplot(DT,aes_string(x="coord",y= "tags",
                                     linetype = "name"))+
            geom_line()+xlab("Genomic Coordinates")+
            ylab(ifelse(normalize,"Normalized Signal","Counts"))
        }else if(tolower(type) == "both"){
          pal = pal[-1]
          out = ggplot(DT,aes_string(x="coord",y= "tags",
                                     linetype = "name",
                                     colour = "type"))+
            geom_line()+xlab("Genomic Coordinates")+
            ylab(ifelse(normalize,"Normalized Signal","Counts"))+
            scale_color_manual(values = pal)
        }else{
          out = ggplot(DT,aes_string(x="coord",y= "tags",
                                     linetype = "name",
                                     colour = "type"))+
            geom_line()+xlab("Genomic Coordinates")+
            ylab(ifelse(normalize,"Normalized Signal","Counts"))+
            scale_color_manual(values = pal)
        }
        out

      })








