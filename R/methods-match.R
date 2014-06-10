
# Methods for match class

## Get methods

#' match1
#'
#' @param match object
#' @return list The match of the reads with "+" strand
#' @docType methods
#' @rdname match-methods
setMethods("match1",
  signature = signature(object = "match"),
  definition = function(object)object@match1
)           

#' match2
#'
#' @param match object
#' @return list The match of the reads with "-" strand
#' @docType methods
#' @rdname match-methods
setMethods("match2",
  signature = signature(object = "match"),
  definition = function(object)object@match2
)


