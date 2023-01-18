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
#' @param format \code{character}. Output file type. See \code{\link[raster]{writeRaster}} for different Options. The default format is "GTiff".
#' @param overwrite \code{logical}. If TRUE, existing files will be overwritten.
#' @return A \code{numeric} raster stack with monthly layers of 
#' aggregated data over specificed time period and area.
#' 
#' @examples
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                     pattern="tas.*\\.nc", full.names=TRUE)
#' summariseRaster(files=files[4], startdate=2001, enddate=2010, var="tas")
#' @export summariseRaster
#' @name summariseRaster
summariseRaster <- function(files, startdate=NA, enddate=NA, ext=NA, var,
                            filename1='', filename2='', format="GTiff", overwrite=FALSE){
  yearmon <- function(x){
    format(as.Date(paste("1", lubridate::month(x, label=T), lubridate::year(x)), "%d %b %Y"), "%b %Y")
  }
  
  if(overwrite==FALSE & file.exists(filename1)){
    avg <- raster::stack(filename1)
  } else{
    if(any(is(files) %in% c("RasterLayer", "RasterStack", "RasterBrick"))){
      data <- files
    } else{
      if(length(files) > 1){
        data <- lapply(files, raster::stack)
        data <- raster::stack(data)
      } else{
        data <- raster::stack(files)
      }
    }
    
    mask <- NA
    if(is(ext, "Extent")){
      # ext is already of class Extent
    } else if(is(ext, "SpatialPolygonsDataFrame")){
      mask <- ext
      ext <- raster::extent(ext)
    } else if(is(ext, "RasterLayer")){
      mask <- ext
      ext <- raster::extent(ext)
    } else if(!anyNA(ext)){
      ext <- raster::extent(ext)
    }
    
    # Crop data by extent
    if(is(ext, "Extent")){
      data <- raster::mask(raster::crop(data, ext), ext)
    }  
    
    # Create list of dates and set dates of raster stack
    dates <- as.Date(gsub("X", "", names(data)), format="%Y.%m.%d")
    data <- raster::setZ(data, dates, 'date')
    
    # Define start date
    if(!is.na(startdate) & any(is(startdate) != "Date")){
      startdate <- as.Date(paste0(startdate, "-01-01"))
    }
    
    # Define end date
    if(!is.na(enddate) & any(is(enddate) != "Date")){
      enddate <-  as.Date(paste0(enddate, "-12-31"))
    }
    
    # Subset dataset by start and enddate
    data <- data[[which(raster::getZ(data) >= startdate & raster::getZ(data) <= enddate)]]
    
    if(var %in% c("hurs", "huss", "tas", "sfcWind")){
      avg <- raster::zApply(data, by=zoo::as.yearmon, fun=mean, name='yearmon', na.rm=TRUE)
    } else if(var == "pr"){
      avg <- raster::zApply(data, by=zoo::as.yearmon, fun=sum, name='yearmon', na.rm=TRUE)
    } else if(var == "tasmax"){
      avg <- raster::zApply(data, by=zoo::as.yearmon, fun=max, name='yearmon', na.rm=TRUE)
    } else if(var == "tasmin"){
      avg <- raster::zApply(data, by=zoo::as.yearmon, fun=min, name='yearmon', na.rm=TRUE)
    }
    avg <- raster::zApply(avg, by=function(x) as.numeric(floor(lubridate::month(as.Date(paste0("01", x), format="%d %b %Y")))), 
                          fun=mean, name='months', na.rm=TRUE)
    names(avg) <- month.name
    
    if(filename1 != ""){
      raster::writeRaster(avg, filename=filename1, format=format, overwrite=overwrite)
    }
    if(filename2 != ""){
      cv <- raster::zApply(data, by=zoo::as.yearmon, fun=cv, name='months', na.rm=TRUE)
      cv <- raster::zApply(cv, by=function(x) as.numeric(floor(lubridate::month(x))), 
                           fun=sum, name='months', na.rm=TRUE)
      raster::writeRaster(cv, filename=filename2, format=format, overwrite=overwrite)
    }
  }; raster::removeTmpFiles(h=0.01)
  return(avg)
}