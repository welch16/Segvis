
# Methods for match class

## Get methods

# @rdname match-methods
# @name match1
# @aliases match
setMethods("match1",
  signature = signature(object = "match"),
  definition = function(object)object@match1
)           

# @rdname match-methods
# @name match2
# @aliases match
setMethods("match2",
  signature = signature(object = "match"),
  definition = function(object)object@match2
)
