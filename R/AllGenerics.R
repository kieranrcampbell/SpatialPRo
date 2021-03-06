
## For SPData

#######################
## Set / Get methods ##
#######################


#' Cell measurement data
#'
#' Returns the n by p matrix of normalised cell measurements for
#' n cells and p channels.
#'
#' @param object The SPData object to use.
#'
#' @name cells
#' @rdname cells-methods
#' @exportMethod cells
#' @docType methods
setGeneric(name = "cells",
           def = function(object) standardGeneric("cells"))

#' Cell measurement data
#'
#' Set the cell measurements for a given sample.
#'
#' @param value The new cell measurements.
#'
#' @name cells<-
#' @rdname cells-methods
#' @exportMethod cells<-
setGeneric(name = "cells<-",
           def = function(object, value) standardGeneric("cells<-"))

#' Raw data matrix
#'
#' Returns the n by p matrix of raw measurements for n cells
#' and p channels.
#'
#' @param object The SPData object to use
#'
#' @name rawData
#' @rdname raw-methods
#' @exportMethod rawData
setGeneric(name = "rawData",
           def = function(object) standardGeneric("rawData"))

#' Number of cells in the sample
#'
#' @param object The SPData object to use
#'
#' @name nCells
#' @rdname nCells-methods
#' @exportMethod nCells
setGeneric(name = "nCells",
           def = function(object) standardGeneric("nCells"))

#' Number of channels in the sample
#'
#' Returns the number of channels, which can be proteins and
#' protein modifications in the given sample.
#'
#' @param object The SPData object to use
#'
#' @name nChannel
#' @rdname nChannel-methods
#' @exportMethod nChannel
setGeneric(name = "nChannel",
           def = function(object) standardGeneric("nChannel"))

#' Channel names used in the sample
#'
#' Returns the names of the channels used in the experiment. This
#' can be the (abbreviated) names of proteins and protein modifications,
#' or whatever else was measured in a given experiment.
#'
#' @param object The SPData object to use
#'
#' @name channels
#' @rdname channels-methods
#' @exportMethod channels
setGeneric(name = "channels",
           def = function(object) standardGeneric("channels"))

#' List of nearest neighbour readouts
#'
#' Returns a list of length nCells(object).
#' The ith entry is an n by m matrix, for cell i having n neighbours
#' each of which have nChannel(object) channels.
#'
#' @param object The SPData object to use
#'
#' @name neighbours
#' @rdname neighbours-methods
#' @exportMethod neighbours
setGeneric(name = "neighbours",
           def = function(object) standardGeneric("neighbours"))

#' List of nearest neighbour readouts
#'
#' Set the list of nearest neighbour readouts.
#'
#' @param value The new list of nearest neighbour readouts
#'
#' @name neighbours<-
#' @rdname neighbours-methods
#' @exportMethod neighbours<-
setGeneric(name = "neighbours<-",
           def = function(object, value) standardGeneric("neighbours<-"))

#' Cell sizes
#'
#' Returns a vector of cell sizes,
#'
#' @param object The SPData object to use
#'
#' @name size
#' @rdname size-methods
#' @exportMethod size
setGeneric(name = "size",
           def = function(object) standardGeneric("size"))

#' Boundary weights
#'
#' Returns a list of boundary weights. The \emph{i}th item
#' will be a vector of length \emph{n} if cell \emph{i} has
#' \emph{n} nearest neighbours.
#'
#' @name weight
#' @rdname weight-methods
#' @exportMethod weight
setGeneric(name = "weight",
           def = function(object) standardGeneric("weight"))

#' Boundary weights
#'
#' Set the boundary weights.
#'
#' @param object The SPData object in which to set the weights
#' @param value A list of length nCells(object) of neighbour weights.
#'
#' @name weight<-
#' @rdname weight-methods
#' @exportMethod weight<-
setGeneric(name = "weight<-",
           def = function(object, value) standardGeneric("weight<-"))

#' Get and set the sample ID
#'
#' Get the sample ID.
#' @param object The SPData object for which to set the ID.
#' @name ID
#' @rdname id-methods
#' @exportMethod ID
setGeneric(name = "ID",
           def = function(object) standardGeneric("ID"))

#' Get and set the sample ID
#'
#' Set the sample ID.
#' 
#' @param value The new sample ID
#'
#' @name ID<-
#' @rdname id-methods
#' @exportMethod ID<-
setGeneric(name = "ID<-",
           def = function(object, value) standardGeneric("ID<-"))

#' Nearest neighbour IDs
#'
#' A list of length \emph{n} for \emph{n} cells, where the \emph{i}th
#' entry is a vector of length \emph{m} that holds the IDs of the nearest
#' neighbour cells of cell \emph{i}, if cell \emph{i} has \emph{m} nearest
#' neighbours.
#'
#' @name neighbourIDs
#' @rdname neighbourid-methods
#' 
#' @param object The \code{SPData} object from which to extract the neighbour IDs.
#' 
#' @exportMethod neighbourIDs
setGeneric(name = "neighbourIDs",
           def = function(object) standardGeneric("neighbourIDs"))

#' 2D Cell Coordinates
#'
#' A matrix of size \emph{n} by 2, for \emph{n} cells. The first
#' column is the \emph{x} coordinate in the tissue, and the
#' second is the \emph{y}.
#'
#' @param object The SPData object from which to retrieve ethe coordinates
#' @rdname xy-methods
#' @exportMethod xy
setGeneric(name = "xy",
           def = function(object) standardGeneric("xy"))

#' 2D Cell Coordinates
#'
#' Set the cell coordinates.
#' 
#' @param value New cell coordinanes. Must be an ncell-by-2 matrix.
#'
#' @name xy<-
#' @rdname xy-methods
#' @exportMethod xy<-
setGeneric(name = "xy<-",
           def = function(object, value) standardGeneric("xy<-"))


#################################
## Class and neighbour methods ##
#################################

#' Classes of each cell
#'
#' If cells within a given tissue have been classified into different
#' types then get/set the classes.
#'
#' @param object The SPData object to use
#'
#' @return A numeric vector of cell classes
#' @rdname cellclass-methods
#' @exportMethod cellClass
setGeneric(name = "cellClass",
           def = function(object) standardGeneric("cellClass"))

#' Classes of each cell
#' 
#' @param value A new vector of cell classes. Must be the same length as the number of cells.
#'
#' @rdname cellclass-methods
#' @name cellClass<-
#' @exportMethod cellClass<-
setGeneric(name = "cellClass<-",
           def = function(object, value) standardGeneric("cellClass<-"))

#' Nearest neighbours of a certain class
#'
#' Filter out nearest neighbours by cell class. If cell \emph{i} has no
#' nearest neighbours of class \code{cell.class} then \code{numeric(0)} is returned.
#'
#' @param object The SPData object to use
#' @param cell.class The cell class to filter out
#'
#' @return A list of length \code{nCell(object)} where all nearest
#' neighbour cells other than of type \code{cell.class} have been removed.
#' @rdname neighbourclass-methods
#' @exportMethod neighbourClass
setGeneric(name = "neighbourClass",
           def = function(object, cell.class)
           standardGeneric("neighbourClass"))







####################
## SPExp generics ##
####################

#' Experiment directory
#'
#' Returns the directory where the experiment files are located.
#'
#' @param object The instance of SPExp to use
#'
#' @rdname getdir-methods
#' @exportMethod getDir
setGeneric(name = "getDir",
           def = function(object) standardGeneric("getDir"))

#' File names
#'
#' Returns a vector of filenames corresponding to the originam experiment files
#'
#' @param object The instance of SPExp to use
#'
#' @rdname files-methods
#' @exportMethod files
setGeneric(name = "files",
           def = function(object) standardGeneric("files"))

#' Underlying list of SPData
#'
#' Accesses the list contained within a \code{SPExp} object, of
#' class \code{SPData}
#'
#' @param object The instance of SPExp to use
#'
#' @rdname splist-methods
#' @exportMethod SPlist
setGeneric(name = "SPlist",
           def = function(object) standardGeneric("SPlist"))

#' Sample IDs
#'
#' Find the IDs of the samples contained within the \code{SPExp}
#'
#' @param object The instance of SPExp to use
#'
#' @rdname ids-methods
#' @exportMethod IDs
setGeneric(name = "IDs",
           def = function(object) standardGeneric("IDs"))

#' Load an experiment
#'
#' Provided the directory and file information exists, this function
#' calls \code{\link{loadCells}} on the files and creates an SPData object.
#'
#' @param object The instance of SPExp to use
#'
#' @return A \code{SPExp} object loaded from \code{getDir(SPExp)}
#'
#' @rdname loadexp-methods
#' @exportMethod loadExp
setGeneric("loadExp", function(object) standardGeneric("loadExp"))

