
# generic methods for various classes

#' @title name
#' @description Generic method for both profile and profileMatrix classes
#' @details This method returns the name of the object
#' @param object Either a profile or profileMatrix object
#' @return character The name of the object
#' @export
#' @docType methods
#' @seealso \code{\link{setName}},\code{\link{profile-class}} and \code{\link{profileMatrix-class}} 
#' @rdname name
setGeneric("name",
  function(object,...)
  standardGeneric("name")
)

#' @title regions
#' @description Generic method for both profile and profileMatrix classes
#' @details This methods returns the regions for which the profile and profileMatrix are calculated as a
#' GRanges object with the regions for which the profiles are calculated
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

#' @title setName
#' @description Generic set name method for both profile and profileMatrix classes
#' @details This methods returns either a profile or profileMatrix object
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

#' @title setRegions
#' @description Generic set regions method for both profile and profileMatrix classes
#' @details This methods returns either a profile or profileMatrix object
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

#' @title file
#' @description Returns the name of the file where the reads are stored
#' @details This method return the value of the file slot
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

#' @title fileFormat
#' @description Returns the format of the file where the reads are stored
#' @details This method return the value of the fileFormat slot
#' @param object A profile object
#' @return character. A character with the file format of the file where reads are stored
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}}
#' @rdname fileFormat
setGeneric("fileFormat",
  function(object,...)
  standardGeneric("fileFormat")         
)

#' @title maxBandwidth
#' @description Returns the maximum bandwidth value
#' @details Return the maxumum bandwidth value
#' @param object A profile object
#' @return numeric value. This indicated the maximum bandwidth available to later smooth the profiles
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{Profile}}
#' @rdname maxBandwidth
setGeneric("maxBandwidth",
  function(object,...)
  standardGeneric("maxBandwidth")           
)

#' @title fragLen
#' @description Returns the fragment length value
#' @details This value indicates the length used to extend the fragment when matching reads with regions
#' @param object A profile object
#' @return numeric value. This indicated the fragment length
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{Profile}}
#' @rdname fragLen
setGeneric("fragLen",
  function(object,...)
  standardGeneric("fragLen")           
)

#' @title remChr
#' @description Returns the chromosomes to be ignored 
#' @details This value indicater the chromosomes that should be ignored
#' @param object A profile object
#' @return character value. Indicates the chromosomes to be ignored
#' @export
#' @docType methods
#' @seealso \code{\link{profile-class}} and \code{\link{Profile}}
#' @rdname remChr
setGeneric("remChr",
  function(object,...)
  standardGeneric("remChr")
)           

#' profileCurve
#' @param profile object
#' @return RleList. A list made of an Rle object for each region
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("profileCurve",
  function(object,...)
  standardGeneric("profileCurve")
)  
    
## Set methods

#' setMaxBandwidth
#'
#' @param profile object
#' @param newMaxBandwidth Numeric value, must be odd and greater or equal than one
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setMaxBandwidth",
  function(object,newMaxBandwidth)
  standardGeneric("setMaxBandwidth")
)           

#' setFragLen
#'
#' @param profile object
#' @param newFragLen Numeric value, must be greater of equal to zero
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setFragLen",
  function(object,newFragLen)
  standardGeneric("setFragLen")
)           

#' setRemChr
#'
#' @param profile object
#' @param newRemChr Character vector, with the chromosomes to be removed
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("setRemChr",
  function(object,newRemChr)
  standardGeneric("setRemChr")
)           

#' loadReads
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("loadReads",
  function(object,mc)
  standardGeneric("loadReads")
)           

#' matchReads
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("matchReads",
  function(object,mc)
  standardGeneric("matchReads")
)           

#' getCoverage
#'
#' @param profile object
#' @param mc numeric, the number of cores used with parallel
#' @return profile object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("getCoverage",
  function(object,mc)
  standardGeneric("getCoverage")           
)

#' buildProfileMat
#'
#' @param profile object
#' @param bw, the bandwidth used to smooth the profiles
#' @param mc, the number of cores used with parallel
#' @return list object
#' @export
#' @docType methods
#' @rdname profile-methods
setGeneric("buildProfileMat",
  function(object,bw,mc)
  standardGeneric("buildProfileMat")
)

# reads generic methods

#' @title reads1
# '@description Get methods for reads class with "+" strand
#' @details  This method return a GRangesList with the reads corresponding to the "+" strand
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

#' @title reads2
#' @description Get methods for reads class with "-" strand
#' @details This method return a GRangesList with the reads corresponding to the "-" strand
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

#' @title setReads1
#' @description Set method of the reads class, for the reads with "+" strand
#' @details This method return a reads object where the reads with "+" strand have been modified by the user
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

#' @title setReads2
#' @description Set method for the reads class, for the reads with "-" strand
#' @details This method return a reads object where the reads with "-" strand have been modified by the user
#' @param reads object
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

#' @title match1
#' @description Get method for match class with "+" strand
#' @details This methods returns a list with the indexes in the "+" strand GRangesList
#' @param match object
#' @return list The match of the reads with "+" strand
#' @export
#' @docType methods
#' @seealso \code{\link{match_reads}},\code{\link{reads1}} and \code{\link{match2}}
#' @rdname match1
setGeneric("match1",
  function(object)
  standardGeneric("match1")
)           

#' @title match2
#' @description Get method for match class with "-" strand
#' @details This methods returns a list with the indexes in the "-" strand GRangesList
#' @param match object
#' @return list The match of the reads with "-" strand
#' @export
#' @docType methods
#' @seealso \code{\link{match_reads}},\code{\link{reads2}} and \code{\link{match1}}
#' @rdname match2
setGeneric("match2",
  function(object,...)
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

            
