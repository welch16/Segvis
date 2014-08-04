
# Methods for match class

## Get methods

# @rdname match-methods
# @name matchF
# @aliases match
setMethod("matchF",
  signature = signature(object = "match"),
  definition = function(object)object@matchF
)           

# @rdname match-methods
# @name matchR
# @aliases match
setMethod("matchR",
  signature = signature(object = "match"),
  definition = function(object)object@matchR
)
