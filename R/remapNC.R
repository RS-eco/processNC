#' Remap NetCDF file
#' 
#' @param gridfile \code{character}. Filepath to the grid.txt file.
#' @param infile \code{character}. Filepath, which should lead to NetCDF files.
#' @param outfile \code{character}. Output filename of the remapped data.
#' @return A NetCDF file containing the croped data.
#' @examples
#' \dontrun{
#' file <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)[1]
#' temp <- tempfile(fileext=".nc")
#' gridfile <- paste0(system.file(package="processNC"), "/extdata/euro-cordex_grid_bav.txt")
#' remapNC(gridfile=gridfile, infile=file, outfile=temp)
#' raster::plot(raster::stack(temp)) 
#' }
#' @export remapNC
#' @name remapNC
remapNC <- function(
  ##title<< Subset data by certain extent
  gridfile ##<< integer vector: lon1, lon2, lat1, lat2
  , infile ##<< character vector: name of the file to crop
  , outfile ##<< character: path to save the results file to. 
)
 {
  ##test input
  #if (system("cdo -V")==0)
  #  stop('cdo not found. Please install it.')
  
  ## supply cdo command
  cdoCmd <- paste(paste('cdo -remapbil', gridfile, sep=","), infile, outfile, sep=' ')
  
  ##run command
  system(cdoCmd)
  cat(paste('Created file ', outfile, '.\n', sep = ''))
  
  ## character string: name of the file created. 
  invisible(outfile)
}
