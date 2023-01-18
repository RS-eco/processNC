#' Crop NetCDF file
#' 
#' @param file \code{character}. Filepath, which should lead to NetCDF files.
#' @param startdate \code{integer}. Start year.
#' @param enddate \code{integer}. End year.
#' @param outfile \code{character}. Output filename of the croped data.
#' @return A NetCDF file containing the croped data.
#' @examples
#' \dontrun{
#' file <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)[4]
#' temp <- tempfile(fileext=".nc")
#' filterNC(file=file, startdate=2005, enddate=2007, outfile=temp)
#' terra::plot(terra::rast(temp))
#' }
#' @export filterNC
#' @name filterNC
filterNC <- function(
  ##title<< Subset data by certain extent
  file ##<< character vector: name of the file to crop
  , startdate
  , enddate 
  , outfile ##<< character: path to save the results file to. 
)
 {
  ##test input
  #if(eval(system("cdo -V")==0))
  #  stop('cdo not found. Please install it.')
  
  ## supply cdo command
  if(length(file)>1){
    cat("Please only provide one file")
  }
  cdoCmd <- paste0('cdo -selyear,', paste0(seq(startdate, enddate), collapse=","),
                   " ", file, ' ', outfile)
  ##run command
  system(cdoCmd)
  cat(paste('Created file ', outfile, '.\n', sep = ''))
  
  ## character string: name of the file created. 
  invisible(outfile)
}
