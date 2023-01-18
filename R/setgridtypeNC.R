#' Set grid type of NetCDF file
#' 
#' @param type \code{character}. Grid type
#' @param infile \code{character}. Filepath, which should lead to NetCDF files.
#' @param outfile \code{character}. Output filename of the remapped data.
#' @return A NetCDF file containing the croped data.
#' @examples
#' \dontrun{
#' file <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)[1]
#' temp <- tempfile(fileext=".nc")
#' setgridtypeNC(type="curvilinear", infile=file, outfile=temp)
#' raster::plot(raster::stack(temp)) 
#' }
#' @export setgridtypeNC
#' @name setgridtypeNC
setgridtypeNC <- function(
  ##title<< Subset data by certain extent
  type="curvilinear" ##<< grid type
  , infile ##<< character vector: name of the file to convert
  , outfile ##<< character: path to save the results file to. 
)
 {
  ##test input
  #if (system("cdo -V")==0)
  #  stop('cdo not found. Please install it.')
  
  ## supply cdo command
  cdoCmd <- paste(paste('cdo -setgridtype', type, sep=","), infile, outfile, sep=' ')
 
  ##run command
  system(cdoCmd)
  cat(paste('Created file ', outfile, '.\n', sep = ''))
  
  ## character string: name of the file created. 
  invisible(outfile)
}
