
##########################
## Utilities for SPData ##
##########################


#' Loads an Xell matlab file into the SPData format
#'
#' This function parses the matlab files, pulling out relevant proteins
#' in the 'D' channel and validates that the correct proteins are present. Matlab
#' files can be found at
#' \url{https://s3.amazonaws.com/supplemental.cytobank.org/report_data/report_113/Figure_5/Figure_5_raw_image_files.zip}
#'
#' @param filename The matlab file
#' @param id The id to give the sample
#' @param control.isotopes The isotopes used for control to exclude from analysis
#'
#' @return An object of class SPData. Note that from the data on cytobank, boundary weights nor cell positions can be found.
#' @export
loadCells <- function(filename,
                      id=-1,
                      control.isotopes = c("Xe131","Cs133","Ir193")) {
  require(R.matlab)
  
  ## loads relevant data from matlab and parses into list
  
  m <- readMat(filename)
  
  n.cells <- dim(m$Xell)[1]
  n.channels <- -1
  
  
  ycolheads <- as.character(m$Xell.list.col)
  yp.id <- grep(")D", ycolheads)
  yp.id <- setdiff(yp.id, grep(paste(control.isotopes, collapse="|"),ycolheads))
  channelNames <- ycolheads[yp.id]
  
  Y <- m$Xell.list
  Y <- Y[,yp.id]
  
  colnames(Y) <- channelNames
  
  ## add minimum to make all values positive
  Y <- preprocess.addmin(Y)
  
  ## fork off 'raw' that this point
  raw <- Y # log(Y)
  
  ## want to LOESS normalise against cell size
  sizes <- as.numeric(m$Xell.size)
  ##Y <- loessNormalise(Y, sizes)
  ##Y <- totalProteinNormalise(Y)
  
  
  ## now on to constructing X, the nearest neighbour matrix
  ## nnids list of nearest neighbour IDs
  nnids <- lapply(m$Xell.nearest, function(xl) {
    if(is.matrix(xl)) {
      return ( xl[,1] )
    } else {
      xl[1]
    }
  })
  
  
  sp <- SPData(channelNames=channelNames,
               readouts=matrix(0), raw=raw, cellNeighbours=list(0),
               size=sizes,id=id, weights=list(0), pos=matrix(0), cellClass=-1,
               nn.ids=nnids)
  return( sp )
}

preprocess.addmin <- function(Y) {
  mu.bg <- -min(Y)
  Y <- Y + mu.bg + 1
}

#' Load an experiment from cytobank in matlab format
#'
#' @param directory The directory containing the experiment files
#' @param files The filenames to load. Default is NULL, in which case all files
#' in the directory are used
#' @return A \code{SPExp} object created from the experiment directory.
#' @export
SPExperimentfromDir <- function(directory, files=NULL) {
  if(is.null(files)) files <- dir(directory)
  
  filesToLoad <- paste(directory,files,sep="/")
  ids <- sapply(filesToLoad, getIDfromTMAname)
  names(ids) <- NULL
  sps <- lapply(1:length(ids), function(i) { loadCells(filesToLoad[i], ids[i])})
  return(SPExperiment(directory, files, sps))
}

getIDfromTMAname <- function(str) {
  splt1 <- strsplit(str, "_ID")[[1]]
  splt2 <- strsplit(splt1[2],"_")[[1]]
  id <- as.numeric(splt2[1])
  id
}


