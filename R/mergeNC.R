#' Merge multiple NetCDF files into one
#' 
#' @param files \code{character}. List of filepaths, which should lead to NetCDF files.
#' @param outfile \code{character}. Output filename of the merged data.
#' @return A NetCDF file containing all of the merged data.
#' @examples
#' \dontrun{
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                     pattern="tas.*\\.nc", full.names=TRUE)
#' temp <- tempfile(fileext=".nc")
#' mergeNC(files=files, outfile=temp)
#' terra::rast(temp) 
#' }
#' @export mergeNC
#' @name mergeNC
mergeNC <- function(
  ##title<< Aggregate data in netCDF files
  files ##<< character vector: names of the files to merge
  , outfile ##<< character: path to save the results files to. 
)
  ##description<<
  ## This function aggregates time periods in netCDF files. Basically it is just a
  ## wrapper around the respective cdo function.
{
  ##test input
  #if (system("cdo -V")==0)
  #  stop('cdo not found. Please install it.')
  
  ## supply cdo command
  cdoCmd <- paste('cdo -cat', paste(files, collapse=" "), outfile, sep=' ')
  
  ##run command
  system(cdoCmd)
  cat(paste('Created file ', outfile, '.\n', sep = ''))
  
  ## character string: name of the file created. 
  invisible(outfile)
}