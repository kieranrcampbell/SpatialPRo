% Generated by roxygen2 (4.0.1): do not edit by hand
\docType{methods}
\name{dim}
\alias{dim}
\alias{dim,SPData-method}
\title{Dimension of underlying matrix representation}
\usage{
\S4method{dim}{SPData}(x)
}
\arguments{
\item{x}{The SPData object to use}
}
\description{
Vector of length 2 that represents the dimensions of the
underlying cell matrix. The first entry is the number of cells
and the second is the number of channels. Equivalent to
dim(cells(x))
}

