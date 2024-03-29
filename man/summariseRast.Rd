% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summariseRast.R
\name{summariseRast}
\alias{summariseRast}
\title{Aggregate values of multiple NetCDF Files}
\usage{
summariseRast(
  files = files,
  startdate = NA,
  enddate = NA,
  ext = NA,
  var,
  filename1 = "",
  filename2 = "",
  overwrite = FALSE
)
}
\arguments{
\item{files}{\code{character}. A filepath or list of filepaths. Filepath should lead to a NetCDF file.}

\item{startdate}{\code{integer}. Start year.}

\item{enddate}{\code{integer}. End year.}

\item{ext}{\code{extent}. If present the NetCDF file is subset by this extent.}

\item{var}{\code{character}. Environmental variable provided.}

\item{filename1}{\code{character}. Output filename of the averaged data. If this argument is not provided, result will not be written to disk.}

\item{filename2}{\code{character}. Output filename of the coefficient of variation. If this argument is not provided, the coefficient of variation will not be written to disk.}

\item{overwrite}{\code{logical}. If TRUE, existing files will be overwritten.}
}
\value{
A \code{numeric} raster stack with monthly layers of 
aggregated data over specificed time period and area.
}
\description{
Aggregate spatio-temporal measurements of one raster stack with daily layers 
for multiple years to a raster stack with monthly layers
}
\examples{
files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
                    pattern="tas.*\\\\.nc", full.names=TRUE)
summariseRast(files=files[4], startdate=2001, enddate=2010, var="tas")
}
