## SPExp - Spatial Proteomics Experiment

#' Class SPExp
#'
#' Class \code{SPExp} defines an entire spatial proteomics experiment.
#' It stores a list of the samples (of type \code{SPData}) as well as
#' the original file names and directories.
#'
#' @name SPExp-class
#' @rdname SPExp-class
#' @aliases SPExp
#'
#' @exportClass SPExp
SPExp <- setClass("SPExp",
                  representation = list(dir = "character",
                      files = "character",
                      spdata = "list"))


#' Displays content of \code{SPExp} object
#' 
#' @aliases show,SPExp-method
#' @export
setMethod("show", "SPExp", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    cat(" Location: ", getDir(object), "\n", sep= "")
    cat( " ", "With ", length(SPlist(object)), " samples \n", sep="")
    cat(" ", IDs(object), "\n", sep=" \n")
    invisible(NULL)
})


#' @rdname getdir-methods
#' @aliases getDir,SPExp-methods
setMethod(f = "getDir",
          signature = "SPExp",
          definition = function(object) object@dir)

#' @rdname files-methods
#' @aliases files,SPExp-methods
setMethod(f = "files",
          signature = "SPExp",
          definition = function(object) object@files)

#' @rdname splist-methods
#' @aliases SPlist,SPExp-methods
setMethod(f = "SPlist",
          signature = "SPExp",
          definition = function(object) object@spdata)

#' @rdname ids-methods
#' @aliases IDs,SPExp-methods
setMethod(f = "IDs",
          signature = "SPExp",
          definition = function(object) sapply(SPlist(object), ID))


#' Initialize SPExp
#' 
#' @param dir The \code{SPExp} directory
#' @param files A list of the original experiment filenames
#' @param spdata A list of \code{SPData} objects
#' 
#' @aliases initialize,SPExp-method
#' @export
setMethod("initialize", signature=signature("SPExp"),
          function(.Object, dir, files, spdata) {
              .Object@dir <- dir
              .Object@files <- files
              .Object@spdata <- spdata
              return(.Object)
          })


#' @rdname loadexp-methods
#' @aliases loadExp,SPExp-methods
setMethod("loadExp", "SPExp",
          function(object) {
              object@spdata <- list()
              N <- length(files(object))
              length(object@spdata) <- N
              for(i in 1:N) {
                  object@spdata[[ i ]] <- loadCells(paste(getDir(object),
                                                          files(object)[i],
                                                          sep=""), id=IDs(object)[i])
              }

              return(object)
          })

#' Subset a \code{SPExp} object.
#'
#' Returns a new \code{SPExp} object subsetted using the first \code{i} instances
#'
#' @param x The \code{SPExp} instance to subset
#' @param i The samples to retain
#' @name [
#' @aliases [,SPExp-method
#' @return A subsetted \code{SPExp} object
#' @rdname spexp-extract-methods
#' @exportMethod [
setMethod("[", "SPExp",
          function(x, i) {
              x@files <- files(x)[i]
              x@spdata <- x@spdata[i]
              return(x)
          })

#' Access a sample
#'
#' Extracts a single sample from a \code{SPExp} object
#'
#' @param x The SPExp instance to use
#' @param i The index of the \code{SPData} object within the \code{SPExp} to extract
#'
#' @return The \code{SPData} object in slot \code{i} of the \code{SPExp} instance
#' @name [[
#' @aliases [[,SPExp-method [[,SPExp,numeric-method
#' @rdname dextract-methods
#' @exportMethod [[
setMethod(f = "[[",
          signature = signature(x="SPExp",i="numeric"),
          function(x,i) {
              return( x@spdata[[i]])
          })


#' Sets a single sample in a \code{SPExp} object
#'
#' @param x The SPExp instance to use
#' @param i The index of the \code{SPData} object within the \code{SPExp} to set
#' @param value An object of class \code{SPData} to replace slot \code{i}
#'
#' @name [[
#' @aliases [[<-,SPExp-method [[<-,SPExp,numeric,missing,SPData-method
#' @rdname spexp-single-extract-methods
setReplaceMethod(f = "[[",
                 signature = signature(x="SPExp",i="numeric",j="missing",value="SPData"),
                 definition = function(x,i, value) {
                     x@spdata[[i]] <- value
                     x
                 })

#' Number of samples
#'
#' @name length
#' @param x \code{SPExp instance to use}
#' @return Number of samples in the \code{SPExp} instance
#' 
#' @aliases length,SPExp-method
#' @export
setMethod(f = "length",
          signature = "SPExp",
          def = function(x) length(x@spdata))

#' Create a new \code{SPExp} object
#' 
#' @param dir Experiment directory
#' @param files Experiment file names
#' @param spdata The list of \code{SPData} objects
#'  
#' @export
SPExperiment <- function(dir, files, spdata) {
  return ( new("SPExp", dir, files, spdata) )
}



#' SPExp object with 5 samples used in report.
#'
#' A dataset of 5 samples of type SPData held in an SPExp object.
#'
#' @docType data
#' @keywords SPExp
#' @name SPE
#' @rdname SPE-data
#' @format An SPExp object with 5 samples
NULL
