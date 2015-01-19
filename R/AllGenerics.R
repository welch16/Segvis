
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

#' @param value The class of this parameter depends on the function, for \code{file<-} and \code{name<-} is a string; for \code{fragLen<-} and  \code{maxBandwidth} is a numeric value and for \code{regions<-} is a GRanges object
#'
#' @return Returns a segvis object with the value indicated by the functions being updated with value
#'
#' @rdname methods-segvis-gs
setGeneric("name<-",
  function(object,value)
  standardGeneric("name<-")
)           

#' @title regions methods for profile and profileMatrix
#' @description Generic set and get methods for both profile and profileMatrix classes that returns the regions for which the profiles are calculated
#' @param object Either a profile or profileMatrix object
#' @return GRanges object with the regions for which the profiles are calculated
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{profileMatrix-class}}
#' @rdname regions
setGeneric("regions",
  function(object)
  standardGeneric("regions")
)

#' @param value GRanges object with the regions for which the profiles are calculated
#' @return profile object with the regions slot replaced by the GRanges object value
#' @export
#' @docType methods
#' @rdname regions
setGeneric("regions<-",
  function(object,value)
  standardGeneric("regions<-")
)

# Generic methods for profile class

##  Get methods

#' @title file method for profile class
#' @description Returns the name of the file where the reads are stored
#' @param object A profile object
#' @return character. A character with the name of the file used to create the reads object
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}
#' @rdname file
setGeneric("file",
  function(object)
  standardGeneric("file")           
)

#' @param value character with the new reads file
#' @return profile object with the file slot replaced by value
#' @export
#' @docType methods
#' @rdname file
setGeneric("file<-",
  function(object,value)
  standardGeneric("file<-")
)

#' @title maxBandwidth method for profile class
#' @description Get and set methods for maxBandwidth parameter of the profile class
#' @param object A profile object
#' @return numeric value. This indicated the maximum bandwidth available to later smooth the profiles
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}
#' @rdname maxBandwidth
setGeneric("maxBandwidth",
  function(object)
  standardGeneric("maxBandwidth")           
)

#' @param newMaxBandwidth Numeric value, must be odd and greater or equal than one
#' @return profile object with the maxBandwidth parameter replaced by value
#' @export
#' @docType methods
#' @rdname maxBandwidth
setGeneric("maxBandwidth<-",
  function(object,value)
  standardGeneric("maxBandwidth<-")
)           

#' @title fragLen get and set methods for profile class
#' @description Returns the fragment length value used to extend the reads when matching them with the regions
#' @param object A profile object
#' @return numeric value. This indicated the fragment length
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}
#' @rdname fragLen
setGeneric("fragLen",
  function(object)
  standardGeneric("fragLen")           
)

#' @param value Numeric, must be greater of equal to zero
#' @return profile object with the fragment length replaced by value
#' @export
#' @docType methods
#' @rdname fragLen
setGeneric("fragLen<-",
  function(object,value)
  standardGeneric("fragLen<-")
)           

#' @title profileCurve methods for profile class
#' @description Returns the coverage for each region considered
#' @param object A profile object
#' @return RleList. A list made of an Rle object for each region
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{getCoverage}}
#' @rdname profileCurve
setGeneric("profileCurve",
  function(object)
  standardGeneric("profileCurve")
)  
    
#' @title loadReads method for profile class
#' @description Load the fragment stored in the file slot of the profile object. The reads are divided by chromosome and by strand
#' @param object profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{readsF}}, \code{\link{readsR}} and \code{\link{reads-class}}
#' @rdname loadReads
#' @examples
#' \dontrun{ loadReads(profile_object,mc=8)}
setGeneric("loadReads",
  function(object,mc)
  standardGeneric("loadReads")
)           

#' @title matchReads methods for profile class
#' @description Match the reads to the extended regions stored in the object slots. 
#' @param object profile object 
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{matchF}}, \code{\link{matchR}} and \code{\link{match-class}}
#' @rdname matchReads
#' @examples
#' \dontrun{ matchReads(profile_object,mc=8)}
setGeneric("matchReads",
  function(object,mc)
  standardGeneric("matchReads")
)           

#' @title getCoverage method for profile class
#' @description Build the coverage using the reads matched by the matchReads method. It returns a Rle object for each region.
#' @param object profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{loadReads}} and \code{\link{matchReads}}
#' @rdname getCoverage
#' @examples
#' \dontrun{ getCoverage(profile_object,mc=8)}
setGeneric("getCoverage",
  function(object,mc)
  standardGeneric("getCoverage")           
)

#' @title buildProfileMatrix method for profile class
#' @description If all the regions have the same width, this method build a profile matrix of dimension nr. regions x width
#' @param object profile object
#' @param bw, the bandwidth used to smooth the profiles, must be and odd number less or equal than the maxBandwidth
#' @param mc, the number of cores used with parallel
#' @return matrix object
#' @export
#' @docType methods
#' @rdname buildProfileMatrix
#' @examples
#' \dontrun{ buildProfileMatrix(profile_object,bw = 151,mc=8)}
setGeneric("buildProfileMatrix",
  function(object,bw,mc)
  standardGeneric("buildProfileMatrix")
)

#' @title findSummit method for profile class
#' @description This method finds the summit of each region whenever is possible, otherwise returns NA
#' @param object profile object
#' @param bw, the bandwidth used to smooth the profiles, must be and odd number less or equal than the maxBandwidth
#' @param mc, the number of cores used with parallel
#' @return Numeric vector
#' @export
#' @docType methods
#' @rdname findSummit
#' @examples
#' \dontrun{ findSummit(profile_object,bw = 151,mc=8)}
setGeneric("findSummit",
  function(object,bw,mc)
  standardGeneric("findSummit")
)

#' @title Create a profileMatrix object
#' @description Constructor for profileMatrix class, it takes both name and region from \code{object}
#' @param object profile object
#' @param bw Numeric used to smooth the profile curves
#' @param mc Numeric, the number of cores used with parallel
#' @return a profileMatrix object
#' @export
#' @docType methods
#' @rdname ProfileMatrix
#' @seealso \code{\link{profileMatrix-class}}
#' @examples
#' \dontrun{ProfileMatrix(profile_object,bw = 151,mc=8)}
setGeneric("ProfileMatrix",
  function(object,bw,mc)
  standardGeneric("ProfileMatrix")
)           

#' @title Generic method readsF for profile and reads classes
#' @description Get method for reads and profile classes that returns the reads from the forward strand
#' @param object reads object
#' @return Returns a GRangesList with the reads of the forward strand
#' @export
#' @seealso \code{\link{reads-class}}
#' @docType methods
#' @rdname readsF
setGeneric("readsF",
  function(object)
  standardGeneric("readsF")
)           

#' @title Generic method readsR for profile and reads classes
#' @description Get method for reads and profile classes that returns the reads from the backward strand
#' @param object reads object
#' @return Returns a GRangesList with the reads of the backward strand
#' @export
#' @seealso \code{\link{reads-class}}
#' @docType methods
#' @rdname readsR
setGeneric("readsR",
  function(object)
  standardGeneric("readsR")           
)

# match generic methods

#' @title matchF method for profile and match classes
#' @description Get method for match class with "+" strand
#' @param match object
#' @return In case of the get method returns the reads of the forward strand
#' @export
#' @docType methods
#' @seealso \code{\link{match-class}}
#' @rdname matchF
setGeneric("matchF",
  function(object)
  standardGeneric("matchF")
)           

#' @title matchR method for profile and match clases
#' @description Get method for match class with "-" strand
#' @param match object
#' @return list The match of the reads with "-" strand
#' @export
#' @docType methods
#' @seealso \code{\link{match-class}}
#' @rdname matchR
setGeneric("matchR",
  function(object)
  standardGeneric("matchR")
)           

# Generic methods for profileMatrix class

## Get methods

#' @title profileMat method for profileMatrix class
#' @description Get and Set methods of the profile matrix for the profileMatrix class 
#' @param object profileMatrix object
#' @return In case of the get method, it returns the profile matrix
#' @export
#' @seealso \code{\link{profileMatrix-class}}
#' @docType methods
#' @rdname profileMat
setGeneric("profileMat",
  function(object)
  standardGeneric("profileMat")           
)        

#' @param value A matrix with a the profile for each region
#' @return In case of the set method, it returns a profileMatrix object where the profile matrix has been replaced by value
#' @export
#' @rdname profileMat
setGeneric("profileMat<-",
  function(object,value)
  standardGeneric("profileMat<-")
)           

#' @title bandwidth method for profileMatrix class
#' @description Get and Set methods for the bandwidth parameter of the profileMatrix class
#' @param object profileMatrix class
#' @return In case of the get method, the value in the bandwidth slot of the object
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#' @rdname bandwidth
setGeneric("bandwidth",
  function(object)
  standardGeneric("bandwidth")           
)           

#' @param value Numeric value with the bandwidth used to smooth the profile in profileMatrix object
#' @return In case of the set method, it returns the profileMatrix object with the bandwidth value replaced by value
#' @export
#' @docType methods
#' @rdname bandwidth
setGeneric("bandwidth<-",
  function(object,value)
  standardGeneric("bandwidth<-")
)           

#' @title normConst method for profileMatrix class
#' @description Get and Set methods for the normConst parameter for profileMatrix class
#' @param object profileMatrix object
#' @return In case of the get method it returns the constant used to normalize the profile
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#' @rdname normConst
setGeneric("normConst",
  function(object)
  standardGeneric("normConst")           
)           

#' @param value Numeric value with the constant used to normalize the profile
#' @return In case of the set method it returns the profileMatrix object with the normalizing constant replced by value
#' @export
#' @docType methods
#' @rdname normConst
setGeneric("normConst<-",
  function(object,value)
  standardGeneric("normConst<-")
)           

#' @title countReads method for profile class
#' @description Counts the number of reads considered  in object
#' @param object profile object
#' @return In case of the get method it returns the constant used to normalize the profile
#' @export
#' @docType methods
#' @seealso \code{\link{profile}}
#' @rdname countReads
setGeneric("countReads",
  function(object)
  standardGeneric("countReads")
)           



#' @title meanProfile method for profileMatrix class
#' @description This method calculates the "average" profile of the matrix respect to their rows. As in the mean function it is possible to set a trim parameter
#' @param object profileMatrix object
#' @param trim Numeric. By default is set to zero. If trim =0, it calculated the row-wise mean, if trim >= 0.5 then is calculates the row-wise median, otherwise it calculate a row-wise trimmed mean.
#' @return Numeric vector
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#' @rdname meanProfile
setGeneric("meanProfile",
  function(object,...)
  standardGeneric("meanProfile")
)           

#' @title normalize.matrix method for profileMatrix class
#' @description This method normalizes the profile matrix by multiplying it times value
#' @param object profileMatrix object
#' @param value A positive number or a missing value. If is missing is going to normalize as if they were 1 million reads in the experiment
#' @return A profileMatrix object where the profileMat slot has been normalized
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#'@rdname normalize.matrix
setGeneric("normalize.matrix",
  function(object,...)
  standardGeneric("normalize.matrix")
)           

#' @title subset.pm method for profileMatrix class
#' @description This method works similarly to the subset of IRanges, GenomicRanges, GenomicAlignments, etc. Althought it doesn't consider the select parameter.
#' @param object profileMatrix object
#' @param condition This is an expression considering the characteristics taht the subset need to satisfy
#' @return Returns a profileMatrix object with the same parameters as object except regions and profileMat which are filtered to satisfy the conditions on condition.
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#' @rdname subset.pm
setGeneric("subset.pm",
  function(object,condition)
  standardGeneric("subset.pm"))

#' @title addColumn method for profileMatrix class
#' @description This method helps to add a new column to the profile matrix object, is works similarly than using the $ operator over \code{regions(object)}
#' @param object profileMatrix object
#' @param name Character with the name of the column to be add
#' @param col Vector of the same lengh as \code{regions(object)}
#' @return Returns a profileMatrix object with the new column added to \code{regions(object)}
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#' @rdname addColumn
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
