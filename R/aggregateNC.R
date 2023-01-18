#' Aggregate one or multiple NetCDF Files over time and/or space
#'
#' @param infile \code{character}. A filepath or list of filepaths. Filepath should lead to a NetCDF file.
#' @param startdate \code{integer}. Start year.
#' @param enddate \code{integer}. End year.
#' @param group_col \code{character}. One of year, yearmon, month
#' @param var \code{character}. Variable which is used, one of tas, tasmin, tasmax and pr.
#' @param outfile \code{character}. Output filename of the averaged data.
#' @return A NetCDF with temporal values averaged over the specified area.
#' @examples
#' \dontrun{
#' file <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)[4]
#' temp <- tempfile(fileext=".nc")
#' aggregateNC(infile=file, outfile=temp, var="tas", startdate="2001", enddate="2010")
#' terra::rast(temp)
#' }
#' @export aggregateNC
#' @name aggregateNC
aggregateNC <- function(
  ##title<< Aggregate data in netCDF files
  infile ##<< character vector: names of the files to aggregate.
  , outfile ##<< character: path to save the results files to. 
  , var ##<< one of tas, tasmax, tasmin, pr
  ,startdate
  ,enddate,
  group_col="year"
  )
  ##description<<
  ## This function aggregates time periods in netCDF files. Basically it is just a
  ## wrapper around the respective cdo function.
{
  
  ##test input
  #if (system("cdo -V")==0)
  #  stop('cdo not found. Please install it.')
  
  # Extract function
  fun <- paste(switch(as.character(var), tasmin='min', tasmax='max', tas='mean', pr='sum'), sep = '')
  
  ## determine cdo command
  if(group_col == "year"){
    cdoCmd <- paste0('cdo -ymonmean -mon', fun, ' -selyear,', 
                     paste0(seq(startdate, enddate), collapse=","), " ", infile, ' ', outfile)
  } else if(group_col == "yearmon"){
    cdoCmd <- paste0('cdo -mon', fun, ' -selyear,', 
                     paste0(seq(startdate, enddate), collapse=","), " ", infile, ' ', outfile)
  } else if(group_col == "month"){
    cdoCmd <- paste0('cdo -year', fun, ' -selyear,', 
                     paste0(seq(startdate, enddate), collapse=","), " ", infile, ' ', outfile)
  }
  
  ##run aggregation
  system(cdoCmd)
  cat(paste('Created file ', outfile, '.\n', sep = ''))
  ##value<<
  ## character string: name of the file created. 
  invisible(outfile)
}