
###################################################
## Contains all plotting routines for SpatialPRo ##
###################################################


#' Boxplots the distribution for each channel
#'
#' @param sp The \code{SPData} object to use
#' @param nrow Number of rows of the boxplots to plot
#' @param ncol Number of columns of boxplots to plot
#'
#' @export
boxplots <- function(sp, nrow, ncol) {
  readouts <- cells(sp)
  read.boundaries <- quantile(readouts, probs=c(0.01, 0.99))
  if(nrow * ncol > nChannel(sp)) {
    print("Not enough channels to boxplot")
    return(NULL)
  } else {
    par(mfrow=c(nrow,ncol), mar=c(1.8,1.8,1.8,1.8))
    for(i in 1:(nrow*ncol)) {
      boxplot(readouts[,i], ylim=read.boundaries,
              main=channels(sp)[i])
    }
  }
}

#' Boxplots the distributions of channels depending on their class
#'
#' @param sp The \code{SPData} object to use
#' @param channel.ids A numeric vector of channel indicies
#' @export
channelPlot <- function(sp, channel.ids) {
  readouts <- cells(sp)
  classes <- as.factor(cellClass(sp))
  readouts <- readouts[,channel.ids]
  colnames(readouts) <- channels(sp)[channel.ids]
  readouts.melted <- melt(readouts)
  plotdf <- data.frame(exprs=readouts.melted$value,
                       channel=readouts.melted$X2,
                       cell.class = rep(classes,length(channel.ids)))
  ggplot(aes(x=channel,y=exprs,fill=cell.class), data=plotdf) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1))
}

#' Correlation plot of channels within the sample
#' 
#' @param sp The \code{SPData} object to di
#' @param cell.class A specific class of cells to pick out for correlation plots
#'
#' @export
channelCorr <- function(sp, cell.class=NULL) {
  Y <- cells(sp)
  if(!is.null(cell.class)) {
    Y <- Y[cellClass(sp) == cell.class,]
  }
  corr.y <- cor(Y)
  corrplot(corr.y)
}

#' Correlation plot of channels with neighbours
#'
#' @param sp The \code{SPData} object to di
#' @param cell.class A specific class of cells to pick out for neighbour correlation plots
#' 
#' @export
neighbourCorr <- function(sp, cell.class=NULL) {
  Y <- cells(sp) ; X <- neighbourMean(sp, TRUE, TRUE)
  if(!is.null(cell.class)) {
    Y <- Y[cellClass(sp) == cell.class,]
    X <- X[cellClass(sp) == cell.class,]
  }
  xy.corr <- cor(X,Y)
  corrplot(xy.corr)
}

