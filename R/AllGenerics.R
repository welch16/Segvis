
# generic methods for various classes

#' @title name method for profile and profileMatrix
#' @description Generic method for both profile and profileMatrix classes that returns the name of the object
#' @param object Either a profile or profileMatrix object
#' @return character The name of the object
#' @export
#' @docType methods
#' @seealso \code{\link{setName}},\code{\link{profile-class}} and \code{\link{profileMatrix-class}} 
#' @rdname name
setGeneric("name",
  function(object)
  standardGeneric("name")
)

#' @title regions methods for profile and profileMatrix
#' @description Generic method for both profile and profileMatrix classes that returns the regions for which the profiles are calculated
#' @param object Either a profile or profileMatrix object
#' @return GRanges object with the regions for which the profiles are calculated
#' @export
#' @docType methods
#' @seealso \code{\link{setRegions}},\code{\link{profile-class}} and \code{\link{profileMatrix-class}}
#' @rdname regions
setGeneric("regions",
  function(object)
  standardGeneric("regions")
)

#' @title setName method for profile and profileMatrix classes
#' @description Generic set name method for both profile and profileMatrix classes
#' @param object Either a profile or profileMatrix object
#' @param newName character
#' @return profile or profileMatrix object
#' @export
#' @docType methods
#' @seealso \code{\link{name}},\code{\link{profile-class}} and \code{\link{profileMatrix-class}} 
#' @rdname setName
setGeneric("setName",
  function(object,newName)
  standardGeneric("setName")           
)

#' @title setRegions methods for profile and profileMatrix classes
#' @description Generic set regions method for both profile and profileMatrix classes
#' @param object Either a profile or profileMatrix object
#' @param newRegions GRanges object with the regions for which the profiles are calculated
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{regions}},\code{\link{profile-class}} and \code{\link{profileMatrix-class}}
#' @rdname setRegions
setGeneric("setRegions",
  function(object,newRegions)
  standardGeneric("setRegions")
)

# Generic methods for profile class

##  Get methods

#' @title file method for profile class
#' @description Returns the name of the file where the reads are stored
#' @param object A profile object
#' @return character. A character with the name of the file used to create the reads object
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{Profile}}
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
#' @seealso \code{\link{profile-class}} and \code{link{Profile}}
#' @rdname fileFormat
setGeneric("fileFormat",
  function(object)
  standardGeneric("fileFormat")         
)

#' @title maxBandwidth method for profile class
#' @description Returns the maximum bandwidth valuee
#' @param object A profile object
#' @return numeric value. This indicated the maximum bandwidth available to later smooth the profiles
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}, \code{\link{Profile}} and \code{\link{setMaxBandwidth}}
#' @rdname maxBandwidth
setGeneric("maxBandwidth",
  function(object)
  standardGeneric("maxBandwidth")           
)

#' @title fragLen method for profile class
#' @description Returns the fragment length value used to extend the reads when matching them with the regions
#' @param object A profile object
#' @return numeric value. This indicated the fragment length
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}, \code{\link{Profile}} and \code{\link{setFragLen}}
#' @rdname fragLen
setGeneric("fragLen",
  function(object)
  standardGeneric("fragLen")           
)

#' @title remChr method for profile class
#' @description Returns the chromosomes to be ignored 
#' @param object A profile object
#' @return character value. Indicates the chromosomes to be ignored
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}, \code{\link{Profile}} and \code{\link{setRemChr}}
#' @rdname remChr
setGeneric("remChr",
  function(object)
  standardGeneric("remChr")
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
    
## Set methods

#' @title setMaxBandwidth method for profile class
#' @description Set method for the maxBandwidth parameter
#' @param object profile object
#' @param newMaxBandwidth Numeric value, must be odd and greater or equal than one
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}, \code{\link{Profile}} and \code{\link{maxBandwidth}}
#' @rdname setMaxBandwidth
setGeneric("setMaxBandwidth",
  function(object,newMaxBandwidth)
  standardGeneric("setMaxBandwidth")
)           

#' @title setFragLen method for profile class
#' @description Set method for the fragLen parameter
#' @param object profile object
#' @param newFragLen Numeric value, must be greater of equal to zero
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}, \code{\link{Profile}} and \code{\link{fragLen}}
#' @rdname setFragLen
setGeneric("setFragLen",
  function(object,newFragLen)
  standardGeneric("setFragLen")
)           

#' @title setRemChr method for profile class
#' @description Set method for the remChr parameter
#' @param object profile object
#' @param newRemChr Character vector, with the chromosomes to be removed
#' @return profile object
#' @export
#' @docType method
#' @seealso \code{\link{profile-class}}, \code{\link{Profile}} and \code{\link{remChr}}
#' @rdname setRemChr
setGeneric("setRemChr",
  function(object,newRemChr)
  standardGeneric("setRemChr")
)           

#' @title loadReads method for profile class
#' @description Load the fragment stored in the file slot of the profile object. The reads are divided by chromosome and by strand
#' @param object profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @seealso \code{\link{reads1}}, \code{\link{reads2}} and \code{\link{reads-class}}
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
#' @seealso \code{\link{match1}}, \code{\link{match2}} and \code{\link{match-class}}
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

#' @title Generic method reads1 for profile and reads classes
# '@description Get methods for reads class with "+" strand
#' @param object reads object
#' @return GRangesList The reads of the "+" strand
#' @export
#' @seealso \code{\link{reads2}}
#' @docType methods
#' @rdname reads1
setGeneric("reads1",
  function(object)
  standardGeneric("reads1")
)           

#' @title Generic method reads2 for profile and reads classes
#' @description Get methods for reads class with "-" strand
#' @param object reads object
#' @return GRangesList The reads of the "-" strand
#' @export
#' @seealso \code{\link{reads1}}
#' @docType methods
#' @rdname reads2
setGeneric("reads2",
  function(object)
  standardGeneric("reads2")           
)

#' @title setReads1 method for reads class
#' @description Set method of the reads class, set the reads in the forward strand
#' @param object reads object
#' @param r1 GRangesList object, this are the new reads to be set in the reads object
#' @return reads object
#' @export
#' @seealso \code{\link{reads1}},\code{\link{reads2}} and \code{\link{setReads2}}
#' @docType methods
#' @rdname setReads1
setGeneric("setReads1",
  function(object,r1)
  standardGeneric("setReads1")
)          

#' @title setReads2 method for reads class
#' @description Set method for the reads class, set the reads in the reverse strand
#' @param object reads object
#' @param r2 GRangesList object, the new reads to set on the reads object
#' @return reads object
#' @export
#' @seealso \code{\link{reads1}},\code{\link{reads2}} and \code{\link{setReads1}}
#' @docType methods
#' @rdname setReads2
setGeneric("setReads2",
  function(object,r2)
  standardGeneric("setReads2")
)          

# match generic methods

#' @title match1 method for profile and match classes
#' @description Get method for match class with "+" strand
#' @param match object
#' @return list The match of the reads with "+" strand
#' @export
#' @docType methods
#' @seealso \code{\link{reads1}} and \code{\link{match2}}
#' @rdname match1
setGeneric("match1",
  function(object)
  standardGeneric("match1")
)           

#' @title match2 method for profile and match clases
#' @description Get method for match class with "-" strand
#' @param match object
#' @return list The match of the reads with "-" strand
#' @export
#' @docType methods
#' @seealso \code{\link{reads2}} and \code{\link{match1}}
#' @rdname match2
setGeneric("match2",
  function(object)
  standardGeneric("match2")
)           

# Generic methods for profileMatrix class

## Get methods

#' profMat
#'
#' @param object. profileMatrix
#' @return matrix object. A matrix with a the profile for each region
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("profMat",
  function(object)
  standardGeneric("profMat")           
)        

#' bandwidth
#'
#' @param object. profileMatrix
#' @return numeric. The bandwidth used to smooth the profile in profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("bandwidth",
  function(object)
  standardGeneric("bandwidth")           
)           

#' normConst
#'
#' @param object. profileMatrix
#' @return numeric. The constant used to normalize the profile
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("normConst",
  function(object)
  standardGeneric("normConst")           
)           

## Set methods

#' setProfMat
#'
#' @param object. profileMatrix
#' @param newProfMat matrix object. A matrix with a the profile for each region
#' @return profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("setProfMat",
  function(object,newProfMatrix)
  standardGeneric("setProfMat")
)           

#' setBandwidth
#'
#' @param object. profileMatrix
#' @param newBandwidth. Numeric value with the bandwidth used to smooth the profile in profileMatrix object
#' @return profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("setBandwidth",
  function(object,newBandwidth)
  standardGeneric("setBandwidth")
)           

#' setNormConst
#'
#' @param object. profileMatrix
#' @param newNormConst Numeric value with the constant used to normalize the profile
#' @return profileMatrix object
#' @export
#' @docType methods
#' @rdname profileMatrix-methods
setGeneric("setNormConst",
  function(object,newNormConst)
  standardGeneric("setNormConst")
)           

            
