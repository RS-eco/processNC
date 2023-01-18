#' Aggregate values of multiple NetCDF Files
#'
#' @description
#' Aggregate spatio-temporal measurements of one raster stack with daily layers 
#' for multiple years to a raster stack with monthly layers
#' 
#' @param files \code{character}. A filepath or list of filepaths. Filepath should lead to a NetCDF file.
#' @param var \code{character}. Environmental variable provided.
#' @param startdate \code{integer}. Start year.
#' @param enddate \code{integer}. End year.
#' @param ext \code{extent}. If present the NetCDF file is subset by this extent.
#' @param filename1 \code{character}. Output filename of the averaged data. If this argument is not provided, result will not be written to disk.
#' @param filename2 \code{character}. Output filename of the coefficient of variation. If this argument is not provided, the coefficient of variation will not be written to disk.
#' @param overwrite \code{logical}. If TRUE, existing files will be overwritten.
#' @return A \code{numeric} raster stack with monthly layers of 
#' aggregated data over specificed time period and area.
#' 
#' @examples
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                     pattern="tas.*\\.nc", full.names=TRUE)
#' summariseRast(files=files[4], startdate=2001, enddate=2010, var="tas")
#' @export summariseRast
#' @name summariseRast
summariseRast <- function(files=files, startdate=NA, enddate=NA, ext=NA, var,
                          filename1='', filename2='', overwrite=FALSE){
  if(overwrite==FALSE & file.exists(filename1)){
    avg <- terra::rast(filename1)
  } else{
    if(any(is(files) %in% c("SpatRaster"))){
      data <- files
    } else{
      if(length(files) > 1){
        data <- lapply(files, terra::rast)
        data <- terra::rast(data)
      } else{
        data <- terra::rast(files)
      }
    }
    
    mask <- NA
    if(is(ext, "Extent")){
      # ext is already of class Extent
    } else if(is(ext, "SpatVector")){
      mask <- ext
      ext <- terra::ext(ext)
    } else if(is(ext, "SpatRaster")){
      mask <- ext
      ext <- terra::ext(ext)
    } else if(!anyNA(ext)){
      ext <- terra::ext(ext)
    }
    
    # Crop data by extent
    if(is(ext, "Extent")){
      data <- terra::mask(terra::crop(data, ext), ext)
    }  
    
    # Create list of dates
    dates <- as.Date(terra::time(data), format="%Y.%m.%d")
    
    # Define start date
    if(!is.na(startdate) & any(is(startdate) != "Date")){
      startdate <- as.Date(paste0(startdate, "-01-01"))
    }
    
    # Define end date
    if(!is.na(enddate) & any(is(enddate) != "Date")){
      enddate <-  as.Date(paste0(enddate, "-12-31"))
    }
    
    # Subset dataset by start and enddate
    data <- data[[which(terra::time(data) >= startdate & terra::time(data) <= enddate)]]
    index_time <- terra::time(data, format="yearmonths")
    index_date <- zoo::as.yearmon(terra::time(data))
    
    if(var %in% c("hurs", "huss", "tas", "sfcWind")){
      avg <- terra::tapp(data, index=index_time, fun=mean, na.rm=TRUE)
    } else if(var == "pr"){
      avg <- terra::tapp(data, index=index_time, fun=sum, na.rm=TRUE)
    } else if(var == "tasmax"){
      avg <- terra::tapp(data, index=index_time, fun=max, na.rm=TRUE)
    } else if(var == "tasmin"){
      avg <- terra::tapp(data, index=index_time, fun=min, na.rm=TRUE)
    }
    terra::time(avg) <- unique(index_date)
    avg <- terra::tapp(avg, index=terra::time(avg, "months"), fun=mean, na.rm=TRUE)
    terra::time(avg) <- c(1:12)
    names(avg) <- month.name
    
    if(filename1 != ""){
      terra::writeCDF(avg, filename=filename1, zname="months", overwrite=overwrite)
    }
    if(filename2 != ""){
      cv <- terra::tapp(data, index=index_time, fun=\(i) raster::cv(i))
      terra::time(cv) <- unique(index_date)
      cv <- terra::tapp(cv, index=terra::time(cv, "months"), fun=sum, na.rm=TRUE)
      terra::time(cv) <- c(1:12)
      names(cv) <- month.name
      terra::writeCDF(cv, filename=filename2, zname="months", overwrite=overwrite)
    }
  }; terra::tmpFiles(remove=T)
  return(avg)
}
