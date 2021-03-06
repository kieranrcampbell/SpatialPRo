% Generated by roxygen2 (4.0.1): do not edit by hand
\name{loadCells}
\alias{loadCells}
\title{Loads an Xell matlab file into the SPData format}
\usage{
loadCells(filename, id = -1, control.isotopes = c("Xe131", "Cs133",
  "Ir193"))
}
\arguments{
\item{filename}{The matlab file}

\item{id}{The id to give the sample}

\item{control.isotopes}{The isotopes used for control to exclude from analysis}
}
\value{
An object of class SPData. Note that from the data on cytobank, boundary weights nor cell positions can be found.
}
\description{
This function parses the matlab files, pulling out relevant proteins
in the 'D' channel and validates that the correct proteins are present. Matlab
files can be found at
\url{https://s3.amazonaws.com/supplemental.cytobank.org/report_data/report_113/Figure_5/Figure_5_raw_image_files.zip}
}

