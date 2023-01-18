#' Convert grid type of NetCDF file
#'
#' @param basedir \code{character}. Path to a file directory containing one or multiple NetCDF files.
#' @return The function overwrites the existing NetCDF files with a new file with the correct grid type, so no new file is returned.
#' @examples
#' \dontrun{
#' extdir <- paste0(system.file(package="processNC"), "/extdata")
#' convertGrid(basedir=extdir)
#' }
#' @export convertGrid
#' @name convertGrid
convertGrid <- function(basedir){
  files <- list.files(basedir, pattern=".nc", full.names=T)
  for(i in 1:length(files)){
    print(files[i])
    # get grid type of input file
    cdoCmd <- paste(paste('cdo -s sinfov', files[i], sep=" "), "| grep points|awk '{print $3}'", sep=' ')
    gridType <- system(cdoCmd, intern=T)
    
    # only fix if grid type is not lonlat
    if(gridType != "lonlat"){
      print("fix grid from $GRIDTYPE to lonlat ...")
      temp <- tempfile(fileext=".nc")
      cdoCmd <- paste(paste0('cdo -s --history -f nc4c -z zip -setgrid,', system.file(package="processNC"), 
                             "/extdata/grid.txt \\"), files[i], temp, "&& mv", temp, file, sep=" ")
      system(cdoCmd)
    } else{
     print(" grid type ok. next...")
    }
  }
}