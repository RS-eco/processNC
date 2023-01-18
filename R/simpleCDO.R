#' Run simple CDO commands
#' 
#' @param command \code{character}. Filepath to the grid.txt file.
#' @param infile \code{character}. Filepath, which should lead to NetCDF files.
#' @return Output
#' @examples
#' \dontrun{
#' file <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)[1]
#' simpleCDO(command="showname", infile=file)
#' simpleCDO(command="showtimestamp", infile=file)
#' simpleCDO(command="sinfo", infile=file)
#' }
#' @export simpleCDO
#' @name simpleCDO
simpleCDO <- function(
  ##title<< Run certain cdo command on data
  command
  , infile ##<< character vector: name of the file
)
 {
  ##test input
  #if (system("cdo -V")==0)
  #  stop('cdo not found. Please install it.')
  
  ## supply cdo command
  cdoCmd <- paste('cdo', command, infile, sep=' ')
  
  ##run command
  system(cdoCmd)
}
