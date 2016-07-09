##' @importFrom viridis viridis
##' @importFrom cba order.optimal
##' @importFrom IRanges overlapsAny
##' @importFrom S4Vectors DataFrame
##' @import ggplot2
##' @import scales
##' @import GenomeInfoDb
##' @import stats
NULL

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
          signature = signature(object = "SegvizData",region = "GRanges"),
          definition = function(object,region,
                                nameFiles = basename(files(object)),
                                type = "aggr",normalize = TRUE,
                                base = 1e6){

            stopifnot(class(region) %in% c("GRanges","SegvizData"))

            stopifnot(is.character(type),
                      tolower(type) %in% c("aggr","fwd","bwd","both","all"))
            stopifnot(is.logical(normalize))

            if(is.numeric(nameFiles))nameFiles = as.character(nameFiles)

            reg = region

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

##' @rdname DT_region-methods
##' @aliases DT_region
##' @docType methods
##' @exportMethod DT_region
setMethod("DT_region",
          signature = signature(object = "SegvizData",region = "numeric"),
          definition = function(object,region,
                                nameFiles = basename(files(object)),
                                type = "aggr",normalize = TRUE,
                                base = 1e6){

            if(is.numeric(region) & !is.integer(region)){
              region = as.integer(region)
            }

            stopifnot(is.character(type),
                      tolower(type) %in% c("aggr","fwd","bwd","both","all"))
            stopifnot(is.integer(region),length(region) == 1,
                      is.logical(normalize))
            stopifnot(1 <= region , region <= length(object))

            if(is.numeric(nameFiles))nameFiles = as.character(nameFiles)

            reg = object[region]

            DT = DT_region(object,reg,nameFiles,type,normalize,base)

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

            tags = NULL
            coord = NULL
            name = NULL

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
            DT = DT[,FUN(tags,...),by = list(coord,name,type)]
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
                geom_line()+xlab("Distance to anchor")+
                ylab("Normalized Signal")
            }else if(tolower(type) == "both"){
              pal = pal[-1]
              out = ggplot(DT,aes_string(x="coord",y= "tags",
                                         linetype = "name",
                                         colour = "type"))+
                geom_line()+xlab("Distance to anchor")+
                ylab("Normalized Signal")+
                scale_color_manual(values = pal)
            }else{

              out = ggplot(DT,aes_string(x="coord",y= "tags",
                                         linetype = "name",
                                         colour = "type"))+
                geom_line()+xlab("Distance to anchor")+
                ylab("Normalized Signal")+
                scale_color_manual(values = pal)
            }
            out

          })

##' @rdname plot_heatmap-methods
##' @aliases plot_heatmap
##' @docType methods
##' @exportMethod plot_heatmap
setMethod("plot_heatmap",
          signature = signature(object = "SegvizData"),
          definition = function(object,which.cluster = 1,
                                dist_method = "euclidean",
                                clust_method = "complete",
                                nameFiles = basename(files(object)),
                                type = "aggr",
                                base = 1e6,
                                mc.cores = getOption("mc.cores",2L),
                                ...){

            coord = NULL
            tags = NULL
            region = NULL

            stopifnot(is.character(dist_method),
                      tolower(dist_method) %in% c("euclidean","maximum",
                                                 "manhattan","canberra",
                                                 "binary","minkowski"))
            if(tolower(dist_method) == "minkowski"){
              stopifnot("p" %in% names(list(...)))
            }

            stopifnot(is.character(clust_method),
                      tolower(clust_method) %in% c("complete","ward.d",
                                                   "ward.d2","single",
                                                   "average","mcquitty",
                                                   "median","centroid"))

            if(!is.integer(which.cluster))which.cluster = as.integer(which.cluster)

            stopifnot(which.cluster >= 1 , which.cluster <= length(files(object)))
            stopifnot(is.character(type),
                      tolower(type) %in% c("aggr","fwd","bwd"))
            stopifnot(length(unique(width(object))) == 1)

            if(is.numeric(nameFiles))nameFiles = as.character(nameFiles)

            aggr_dt = NULL
            fwd_dt = NULL
            bwd_dt = NULL
            .x = NULL

            width = unique(width(object))
            ll = (width - 1) / 2

            if(tolower(type) == "aggr"){
              aggr_dt = mapply(.dt_profile,covers(object),
                               nreads(object),nameFiles,
                               MoreArgs = list(regions = object,
                                               st = "aggr",
                                               base = base,
                                               len = ll,
                                               mc.cores = mc.cores),
                               SIMPLIFY = FALSE)
            }

            if(tolower(type) == "fwd"){
              fwd_dt = mapply(.dt_profile,fwd_covers(object),
                              nreads(object),nameFiles,
                              MoreArgs = list(regions = object,
                                              st = "fwd",
                                              base = base,
                                              len = ll,
                                              mc.cores = mc.cores),
                              SIMPLIFY = FALSE)

            }
            if(tolower(type) == "bwd"){

              bwd_dt = mapply(.dt_profile,bwd_covers(object),
                              nreads(object),nameFiles,
                              MoreArgs = list(regions = object,
                                              st = "bwd",
                                              base = base,
                                              len = ll,
                                              mc.cores = mc.cores),
                              SIMPLIFY = FALSE)

            }

            dt_list = c(fwd_dt,bwd_dt,aggr_dt)

            to_clust = dt_list[[which.cluster]][,list(coord,tags,region)]
            to_clust = dcast.data.table(data = to_clust,
                                        formula = region ~ coord,
                                        value.var = "tags",
                                        fun.aggregate = mean)
            if(tolower(dist_method ) == "minkowski"){
              dmat = dist(as.matrix(to_clust[,-1,with = FALSE]),
                          method = "minkowski",p = list(...)$p)
            }else{
              dmat = dist(as.matrix(to_clust[,-1,with = FALSE]),
                        method = dist_method)
            }
            clust = hclust(dmat,method = clust_method)

            ord = cba::order.optimal(dmat,clust$merge)

            lev = unique(dt_list[[which.cluster]][,(region)])[ord$order]

            DT = rbindlist(dt_list)
            DT = DT[,region := factor(region,levels = lev)]

            p = ggplot(DT,aes_string(x = "coord",y = "region",fill = "tags"))+
              geom_raster()+
              scale_fill_gradientn(name = "Normalized Signal",
                    colors = viridis(50),trans = 'log10',
                    labels = trans_format('log10',math_format(10^.x)))+
              facet_grid( . ~ name)+
              geom_vline(xintercept = 0,linetype = 2,colour = "lightgrey")+
              theme(legend.position = "top",
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.key.width = unit(2,"lines"))+
              xlab("Distance to anchor")+ylab("")
            p
          })

##' @rdname plot_region-methods
##' @aliases plot_region
##' @docType methods
##' @exportMethod plot_region
setMethod("plot_region",
          signature = signature(object = "SegvizData",region = "GRanges"),
          definition = function(object,region,
                                nameFiles = basename(files(object)),
                                type = "aggr",normalize = TRUE,
                                base = 1e6){

            DT = DT_region(object,region,nameFiles = nameFiles,
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


##' @rdname plot_region-methods
##' @aliases plot_region
##' @docType methods
##' @exportMethod plot_region
setMethod("plot_region",
      signature = signature(object = "SegvizData",region = "numeric"),
      definition = function(object,region,
                            nameFiles = basename(files(object)),
                            type = "aggr",normalize = TRUE,
                            base = 1e6){

        plot_region(object,region = object[region],nameFiles,type,
                    normalize,base)
      })








