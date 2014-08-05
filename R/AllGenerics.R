
# generic methods for various classes

#' @title name get and set methods for profile and profileMatrix
#' @description Generic method for both profile and profileMatrix classes that returns the name of the object
#' @param object Either a profile or profileMatrix object
#' @return Returns the name of the object
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{profileMatrix-class}}
#' @rdname name
setGeneric("name",
  function(object)
  standardGeneric("name")
)

#' @param value character 
#' @return Returns and profile (or profileMatrix) object with the name parameter replaced by value
#' @rdname name
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

#' @param newRegions GRanges object with the regions for which the profiles are calculated
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

#' @title fileFormat method for profile class
#' @description Returns the format of the file where the reads are stored
#' @param object A profile object
#' @return character. A character with the file format of the file where reads are stored
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}
#' @rdname fileFormat
setGeneric("fileFormat",
  function(object)
  standardGeneric("fileFormat")         
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

#' @title remChr get and set methods for profile class
#' @description Returns the chromosomes to be ignored 
#' @param object A profile object
#' @return character value. Indicates the chromosomes to be ignored
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}
#' @rdname remChr
setGeneric("remChr",
  function(object)
  standardGeneric("remChr")
)           

#' @param value Character vector, with the chromosomes to be removed
#' @return profile object with the remChr parameter replaced by value
#' @export
#' @docType method
#' @rdname remChr
setGeneric("remChr<-",
  function(object,value)
  standardGeneric("remChr<-")
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

#' @title subset.pm method for profileMatrix class
#' @description This method works similarly to the subset of IRanges, GenomicRanges, GenomicAlignments, etc. Althought it doesn't consider the select parameter.
#' @param object profileMatrix object
#' @param subset This is an expression considering the characteristics taht the subset need to satisfy
#' @return a profileMatrix object with the same parameters as object except regions and profileMat which are filtered to satisfy the conditions on subset.
#' @export
#' @docType methods
#' @seealso \code{\link{profileMatrix-class}}
#' @rdname subset.pm
setGeneric("subset.pm",
  function(object,subset)
  standardGeneric("subset.pm"))


