
# Neighbour methods -------------------------------------------------------



#' Average over nearest neighbours
#'
#' Given the nearest neighbour measurements, average over the
#' nearest neighbours of each cell, with optional boundary
#' weights.
#'
#' @param object The SPData object to use
#' @param useWeights If TRUE then the mean is weighted using the
#' relative boundary weights. If FALSE then the normal mean is taken
#' over the nearest neighbours.
#' @param normalise If TRUE each channel is centre-scaled to mean 0
#' and standard deviation 1.
#'
#' @return A matrix of dimension \emph{n} by \emph{p} for
#' \emph{n} cells and \emph{p} nearest neighbours.
#' @rdname neighbourmean-methods
#' @export
neighbourMean  <-  function(object, useWeights = FALSE, normalise = TRUE) {
  ## average over nearest neighbours then means
  X <- neighbours(object)
  
  weights <- NULL
  if(useWeights) weights <- weight(object)
  
  X <- lapply(1:length(X), function(i) {
    nn <- X[[i]]
    
    if(is.matrix(nn)) {
      if(useWeights) {
        w <- weights[[i]]
        total.boundary <- sum(w)
        nn <- nn * w / total.boundary ## IMPORTANT: matrix * vector multiplication is by column
        return( colSums(nn) )
      } else {
        return( colMeans(nn) )
      }
    }
    else {
      return( nn )
    }
  })
  
  X <- matrix(unlist(X), nrow=nCells(object), ncol=nChannel(object), byrow=TRUE)

  if(normalise) X <- apply(X, 2, function(x) (x - mean(x))/sd(x))
  
  colnames(X) <- channels(object)
  X
}

#' Selects only particular channels to be returned as nearest neighbours
#'
#' @param NN The list of matrices of neighbour data (e.g. from \code{neighbours(sp)} )
#' @param channel.ids A channel list to select out
#'
#' @export
neighbourChannel <- function(NN, channel.ids) {
  ## the trick with this function is to make sure each element in the nearest neighbour
  ## list arrives as a matrix and leaves as a matrix. R will reduce an
  
  nn <- lapply(NN, function(Xi) {
    if(!is.matrix(Xi)) {
      if(is.na(Xi)) { 
        return( NA ) 
      } else {
        stop("Dimensionality lost in neighbour subsetting")        
      }
    } else {
      Xi <- Xi[,channel.ids, drop = FALSE]
      return( Xi )
    }
  })
  return(nn)
}

#' Generate nearest neighbour data
#' 
#' Some processing may be performed on the cell-by-channel matrix and the nearest neighbour
#' data needs updated. Note that the class architecture means that when the neighbours of a cell
#' are requested, it is not generated fresh from the cell matrix but retreived from a separate copy 
#' specifically for neighbour data. Therefore, if any transformations are applied to the neighbour
#' data, you will need to regenerate them using this function.
#' 
#' @param sp The \code{SPData} object to use
#' @return A \code{SPData} object with the neighbour data regenerated
#' 
#' @export
generateNeighbourfromID <- function(sp) {
  Y <- cells(sp)
  X <- lapply(neighbourIDs(sp), function(id) Y[id, ,drop=FALSE])
  neighbours(sp) <- X
  return ( sp )
}

#' Cells on a class boundary
#'
#' Find the cells that lie on the boundary between two classes
#' (currently only implemented for 2 classes)
#'
#' @param sp The SPData object to use
#'
#' @return A vector of cell identifiers relating to those
#' that lie along the boundary
#'
#' @export
findBoundary <- function(sp) {
  classes <- cellClass(sp)
  cl1 <- which(classes == 1) ; cl2 <- which(classes == 2)
  nn.ids <- neighbourIDs(sp)
  
  cell1neighbours <- sapply(cl1, function(cellid) {
    nn.id <- nn.ids[[cellid]]
    any(nn.id %in% cl2)
  })
  
  cell2neighbours <- sapply(cl2, function(cellid) {
    nn.id <- nn.ids[[cellid]]
    any(nn.id %in% cl1)
  })
  boundary <- sort(c(cl1[cell1neighbours], cl2[cell2neighbours]))
  return( boundary )
}

checkNeighbours <- function(sp) {
  null <- sapply(neighbours(sp), function(nn) {
    if(!is.matrix(nn)) {
      print(nn)
      stop("nn not matrix")
    }
  })
}
