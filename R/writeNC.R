#' Save data to a NetCDF file
#'
#' @param data \code{list}. data.frame or list of data.frames according to the length of the variables.
#' @param res \code{numeric}. Resolution of the output file.
#' @param ext \code{extent}. If present the NetCDF file is subset by this extent.
#' @param date \code{numeric}. Date(s) of input data which should be written to file.
#' @param vars \code{character}. Variables included in output file.
#' @param varunit \code{character}. Unit of variable output.
#' @param contact \code{character}. Contact information of author.
#' @param institute \code{character}. Institute affiliation of author.
#' @param filename \code{character}. Name of output file.
#' @return A raster stack with monthly layers of aggregated data over the specified time period and area.
#' @examples
#' \dontrun{
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                    pattern="tas.*\\.nc", full.names=TRUE)
#' data <- as.data.frame(raster::rasterToPoints(raster::stack(files)))
#' outfile <- tempfile(fileext=".nc4")
#' writeNC(data=data, res=0.5, ext=c(-180,180,-90,90), 
#'         date=seq(as.Date("1979-01-01"), as.Date("2016-12-31"), by="day"), 
#'         vars="tas", varunit="in degree Celsius", contact="RS-eco <rs-eco@posteo.de>", 
#'         institute="Github", filename=outfile)
#' raster::stack(outfile)
#' }
#' @export writeNC
#' @name writeNC
writeNC <- function(data, res, ext, date, vars, 
                    varunit, contact, institute,
                    filename){
  
  # Define the dimensions
  dimX = ncdf4::ncdim_def(name="lon", units="degrees", vals=seq(ext[1]+(res/2), ext[2]-(res/2), length=(ext[2]-ext[1])/res))
  dimY = ncdf4::ncdim_def(name="lat", units="degrees", vals=seq(ext[4]-(res/2), ext[3]+(res/2), length=(ext[4]-ext[3])/res))
  dimT = ncdf4::ncdim_def(name="time", units="days since 1661-1-1 00:00:00", 
                          vals=as.numeric(date-1661), calendar="standard")
  
  # Define data for NetCDF file
  vard <- lapply(vars, function(var){ncdf4::ncvar_def(var, varunit, list(dimX,dimY,dimT), 1.e+20, prec="double", compression=9)})
  
  # Create the NetCDF file
  nc <- ncdf4::nc_create(filename, vard)
  ncdf4::ncatt_put(nc, varid=0, attname="contact", attval=contact)
  ncdf4::ncatt_put(nc, varid=0, attname="institution", attval=institute)
  
  lapply(1:length(vars), function(x){
    var <- vars[x]
    
    ## Select all data for one varible
    if(class(data) == list){
      data_sub <- data[[x]]
    } else{
      data_sub <- data
    }
    
    # Expand dataframe with NAs
    df_spat <- expand.grid(x=seq(ext[1]+(res/2), ext[2]-(res/2), length=(ext[2]-ext[1])/res),
                           y=seq(ext[4]-(res/2), ext[3]+(res/2), length=(ext[4]-ext[3])/res))
    data <- dplyr::left_join(df_spat, data_sub); rm(df_spat)
    
    # Turn data into array
    data_sub <- dplyr::select(data_sub, -"x", -"y")
    colnames(data_sub) <- date
    data <- array(unlist(data_sub), dim=c((ext[2]-ext[1])/res, (ext[4]-ext[3])/res, ncol(data_sub)), 
                  dimnames=list(NULL, NULL, colnames(data_sub)))
    
    # Write data to the NetCDF file
    ncdf4::ncvar_put(nc, vard[[x]], data, start=c(1,1,1), count=c(-1,-1,-1))
  })
  
  # Close your new file to finish writing
  ncdf4::nc_close(nc); rm(nc)
}
