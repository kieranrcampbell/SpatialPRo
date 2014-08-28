################################################
## SpatialPRo R class                         ##
## kieran.campbell@dpag.ox.ac.uk              ##
## Data container for spatial proteomics data ##
################################################

#' Class SPdata
#'
#' Class \code{SPData} defines a spatial proteomics sample. It
#' stores a sample (e.g. a tissue microarray) of a given
#' tissue as single cell proteomics data. It contains
#' cell measurements, size and classification.
#'
#' @name SPData-class
#' @rdname spdata-class
#' @aliases SPData
#'
#' @exportClass SPData
SPData <- setClass(Class = "SPData",
                   representation = list(channelNames = "character",
                       readouts = "matrix", ## +min
                       raw = "matrix", ## + min + log
                       cellNeighbours = "list",
                       nn.ids = "list",
                       size = "numeric",
                       id = "numeric",
                       weights = "list",
                       pos = "matrix", # nCell by 2 matrix of cell locations
                       cellClass = "numeric"))
#' Some text
#' @name cells
#' @rdname cells-methods
#' @aliases cells,SPData-method
setMethod(f = "cells",
          signature = signature(object="SPData"),
          definition = function(object) object@readouts )

#' @rdname cells-methods
#' @name cells
#' @aliases cells<-,SPData,matrix-method cells<-,SPData-method
setReplaceMethod("cells", signature(object="SPData",value="matrix"),
                 function(object, value) {
                     object@readouts <- value
                     validObject(object)
                     return(object)
                 })

#' @rdname raw-methods
#' @aliases rawData,SPData-method
setMethod(f = "rawData",
          signature = "SPData",
          definition = function(object) object@raw)

#' @rdname nCells-methods
#' @aliases nCells,SPData-method
setMethod(f = "nCells",
          signature = "SPData",
          definition = function(object) dim(object@raw)[1] )

#' @rdname nChannel-methods
#' @aliases nChannel,SPData-method
setMethod(f = "nChannel",
          signature = "SPData",
          definition = function(object) dim(object@readouts)[2] )

#' @rdname channels-methods
#' @aliases channels,SPData-method
setMethod(f = "channels",
          signature = "SPData",
          definition = function(object) object@channelNames )

#' @rdname neighbours-methods
#' @aliases neighbours,SPData-method
setMethod(f = "neighbours",
          signature = "SPData",
          definition = function(object) {
    object@cellNeighbours
})

#' @rdname neighbours-methods
#' @name neighbours
#' @aliases neighbours<-,SPData-method neighbours<-,SPData,list-method
setReplaceMethod("neighbours", signature = signature(object="SPData",value="list"),
                 function(object, value) {
                     object@cellNeighbours <- value
                     validObject(object)
                     return(object)
                 })


#' @rdname size-methods
#' @aliases size,SPData-method
setMethod(f = "size",
          signature = "SPData",
          definition = function(object) object@size)

#' @rdname weight-methods
#' @aliases weight,SPData-method
setMethod(f = "weight",
          signature = "SPData",
          def = function(object) object@weights)

#' @rdname weight-methods
#' @name weight<-
#' @aliases weight<-,SPData-method weight<-,SPData,list-method
setReplaceMethod(f = "weight",
                 signature=signature(object="SPData",value="list"),
                 function(object, value) {
                     object@weights <- value
                     return(object)
                 })


#' Dimension of underlying matrix representation
#'
#' Vector of length 2 that represents the dimensions of the
#' underlying cell matrix. The first entry is the number of cells
#' and the second is the number of channels. Equivalent to
#' dim(cells(x))
#'
#' @param x The SPData object to use
#' @name dim
#' @rdname dim-methods
#' @aliases dim,SPData-method
#' @exportMethod dim
setMethod(f = "dim",
          signature = "SPData",
          def  = function(x) c(nCells(x), nChannel(x)))

#' @rdname id-methods
#' @aliases ID,SPData-method
setMethod(f = "ID",
          signature = "SPData",
          def  = function(object) object@id)

#' @name ID<-
#' @rdname id-methods
#' @aliases ID<-,SPData-method ID<-,SPData,numeric-method
setReplaceMethod(f = "ID",
                 signature = signature(object="SPData", value="numeric"),
                 function(object, value) {
                     object@id <- value
                     return(object)
                 })

#' @rdname neighbourid-methods
#' @aliases neighbourIDs,SPData-method
setMethod(f = "neighbourIDs",
          signature = "SPData",
          def = function(object) object@nn.ids)

#' @rdname xy-methods
#' @aliases xy,SPData-method
setMethod(f = "xy",
          signature = "SPData",
          def = function(object) object@pos)


#' @rdname xy-methods
#' @aliases xy<-,SPData,matrix-method xy<-,SPData-method
#' @name xy<-
#' @exportMethod xy<-
setReplaceMethod(f = "xy",
                 signature=signature(object="SPData",value="matrix"),
                 function(object, value) object@pos <- xy)



#' Default show call
#' @export
setMethod("show", "SPData", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    if(ID(object) > -1) cat(" Sample ID: ", ID(object), "\n", sep="")
    cat(" ", nCells(object), " cells with ",
        nChannel(object), " channel(s)\n", sep="")
    invisible(NULL)

})


setValidity("SPData", function(object) {
    msg <- NULL
    valid <- TRUE
    if(length(neighbours(object)) > 1) {
        if(nCells(object) != length(neighbours(object))) {
            valid <- FALSE
            msg <- c(msg, "Nearest neighbour data not available for all cells")
        }
    }

    if(!is.null(xy(object)) && length(xy(object)) > 1) {
        if(nCells(object) != nrow(xy(object))) {
            valid <- FALSE
            msg <- c(msg, "Length mismatch between location info and number of cells")
        }
    }

    if(length(weight(object)) > 1 ) { ## okay for object not to have weights
        if(length(weight(object)) != nCells(object)) {
            valid <- FALSE
            msg <- c(msg, "Number of weights and number of cells differ")
        }
    }

    if(nCells(object) != nrow(rawData(object))) {
        valid <- FALSE
        ##print(nCells(object))
        ##print(nrow(cells(object)))
        msg <- c(msg, "Number of cells must be equal to number of rows in cell by protein matrix")
    }

    if(nChannel(object) != dim(cells(object))[2]) {
        valid <- FALSE
        msg <- c(msg, "Number of proteins must be equal to number of columns in cell by protein matrix")
    }
    
    if(length(cellClass(object)) > 0 && length(cellClass(object)) != nCells(object)) {
      valid  <- FALSE
      msg  <- c(msg, "Length of cell class vector doesn't match number of cells")
    }

    if(valid) TRUE else msg

})


#' Subset an SPData set
#'
#' Select SPData[i,j] for cells \code{i} and channels \code{j}.
#' Note that this does not subset out nearest neighbours also.
#' 
#' @name [
#'
#' @param i Cells to subset
#' @param j Channels to subset
#'
#' @return An SPData object reduced to cells \code{i} and channels \code{j}
#'
#' @aliases [,SPData-method [,SPData,ANY,ANY,ANY-method
#' 
#' @rdname extract-methods
#' @docType methods
#' @exportMethod [
#' @examples 
#' \dontrun{
#' ## subset to cells 1,3,5 and channels 8 to 10:
#' i <- c(1,3,5)
#' j <- 8:10
#' sp.reduced <- sp[i,j]
#' }
#' 
setMethod("[", signature(x="SPData",i="ANY",j="ANY"), function(x, i, j) {
    if(missing(j)) j <- 1:nChannel(x)
    if(missing(i)) i <- 1:nCells(x)
    
    if(is.logical(i)) i <- which(i)
    if(is.logical(j)) j <- which(j)

    .n.proteins <- length(j)
    .channelNames <- channels(x)[j]
    .id <- ID(x)
    .weight <- weight(x)[i]
    .nnid <- neighbourIDs(x)[i]
    .cell.class <- cellClass(x)[i]
    .pos <- xy(x)[i,,drop=FALSE]

    ## if dealing with single cell, .pos will become vector:
    if(!is.matrix(.pos)) .pos <- t(as.matrix(.pos))
    
    .Y <- cells(x)[i,j, drop=FALSE]
    .raw <- rawData(x)[i,j, drop=FALSE]

    .X <- neighbourChannel(neighbours(x), j)

    .X <- .X[i]

    .size <- size(x)[i]

    SPData(channelNames = .channelNames,
           readouts = .Y,
           cellNeighbours = .X,
           size = .size, id=.id,
           weights = .weight, nn.ids = .nnid,
           raw = .raw, cellClass = .cell.class, pos=.pos)
})

###########################################
## Methods for neighbours and cell class ##
###########################################


#' @rdname cellclass-methods
#' @aliases cellClass,SPData-method
setMethod("cellClass", signature="SPData", function(object) object@cellClass)

#' @name cellClass<-
#' @rdname cellclass-methods
#' @aliases cellClass<-,SPData-method cellClass<-,SPData,numeric-method
setReplaceMethod("cellClass", signature=signature(object="SPData", value="numeric"),
                 function(object, value) {
                   if(length(value) != nCells(object)) stop("Length of class vector different to number of cells in SPData object")                   
                   object@cellClass <- value
                   return(object)
                 })

#' @rdname neighbourclass-methods
#' @aliases neighbourClass,SPData-method
setMethod("neighbourClass", signature("SPData","numeric"),
          function(object, cell.class) {
              if(!(cell.class %in% cellClass(object))) stop("Cell class not present in tissue")
              X <- neighbours(object)
              nn.ids <- neighbourIDs(object)
              cell.select <- which(cellClass(object) == cell.class)

              nn <- lapply(1:length(X), function(i) {
                  Xi <- X[[i]] ; ids <- nn.ids[[i]]

                  lvec <- ids %in% cell.select
                  if(is.matrix(Xi)) {
                    if(any(lvec)) return(Xi[lvec,,drop=FALSE]) else return( NA )
                  } else {
                    if(is.na(Xi)) return( NA ) else stop("Nearest neighbour matrix must be matrix or NA")  
                  }
              })
              nn
          })











