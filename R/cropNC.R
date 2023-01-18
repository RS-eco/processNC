#' Crop NetCDF file
#' 
#' @param file \code{character}. Filepath, which should lead to NetCDF files.
#' @param ext \code{integer}. Extent object.
#' @param outfile \code{character}. Output filename of the croped data.
#' @return A NetCDF file containing the croped data.
#' @examples
#' \dontrun{
#' file <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)[1]
#' temp <- tempfile(fileext=".nc")
#' cropNC(file=file, ext=c(9, 13, 49, 51), outfile=temp)
#' terra::plot(terra::rast(temp)) 
#' }
#' @export cropNC
#' @name cropNC
cropNC <- function(
  ##title<< Subset data by certain extent
  file ##<< character vector: name of the file to crop
  , ext ##<< integer vector: lon1, lon2, lat1, lat2
  , outfile ##<< character: path to save the results file to. 
)
 {
  ##test input
  #if(eval(system("cdo -V")==0))
  #  stop('cdo not found. Please install it.')
  
  ## supply cdo command
  cdoCmd <- paste(paste('cdo -selindexbox', ext[1], ext[2], ext[3], ext[4], sep=","), file, outfile, sep=' ')
  
  ##run command
  system(cdoCmd)
  cat(paste('Created file ', outfile, '.\n', sep = ''))
  
  ## character string: name of the file created. 
  invisible(outfile)
}
