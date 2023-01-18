#' Subset NetCDF files by x,y,time or variable
#'
#' @description
#' Subset a NetCDF file by a spatial extent (x,y), time or variable and 
#' output it as a raster stack 
#'  
#' @param files \code{character}. File path.
#' @param ext \code{integer}. Extent object.
#' @param startdate \code{integer}. Start year.
#' @param enddate \code{integer}. End year.
#' @param var \code{character}. Variable that should be extracted.
#' @param filename \code{character}. Output filename.
#' @param format \code{character}. Output file type. See \code{\link[raster]{writeRaster}} for different Options. The default format is "raster".
#' @param overwrite \code{logical}. Should file be overwritten or not.
#' @return A \code{numeric} raster stack with the cropped NetCDF data.
#' @examples
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                     pattern="tas.*\\.nc", full.names=TRUE)
#' subsetNC(files, ext=c(8.5, 14, 47, 51), startdate=1990, enddate=1999)
#' @export subsetNC
#' @name subsetNC
subsetNC <- function(files, startdate=NA, enddate=NA, ext="",
                     var="", filename="", format="raster", overwrite=FALSE){
  
  if(overwrite==FALSE & filename != "" & file.exists(filename)){
    r <- terra::rast(filename)
  } else{
    mask <- NA
    if(any(is(ext) != "Extent")){
      if(is(ext, "SpatVector")){
        mask <- ext
        ext <- terra::ext(ext)
      } else if(is(ext, "SpatRaster")){
        ext <- terra::ext(ext)
      } else if(any(ext != "")){
        ext <- terra::ext(ext)
      } else{
        # Obtain size and other specifications of file
        nc <- ncdf4::nc_open(files[1])
        
        # Obtain coordinates
        lon <- ncdf4::ncvar_get(nc, nc$dim$lon, start=1)
        lat <- ncdf4::ncvar_get(nc, nc$dim$lat, start=1)
        
        # Create extent object
        ext <- c(min(lon), max(lon), min(lat), max(lat))
        ext <- terra::ext(ext)
        
        # Close NC file again
        ncdf4::nc_close(nc)
      }
    }
    
    # Convert start and endyear to dates
    if(!is.na(startdate) & any(is(startdate) != "Date")){startdate <- as.Date(paste0(startdate, "-01-01"))}
    if(!is.na(enddate) & any(is(enddate) != "Date")){enddate <- as.Date(paste0(enddate, "-01-01"))}
    
    # Obtain size and other specifications of file
    data <- lapply(files, FUN=function(a){
      nc <- ncdf4::nc_open(a)
      
      #Set count and start to	1,1,1,...,1
      count <- rep(1,nc$ndims)
      start <- rep(1,nc$ndims)
      
      # Get coordinates of file
      x <- ncdf4::ncvar_get(nc, nc$dim$lon)
      y <- ncdf4::ncvar_get(nc, nc$dim$lat)
      
      # Specify varsize and number of dimensions
      varsize <- nc$var[[nc$nvars]]$varsize
      ndims   <- nc$var[[nc$nvars]]$ndims
      
      # Get start and end date of file
      timeref <- as.Date(strsplit(nc$dim[[nc$ndims]]$units, " ")[[1]][3]) 
      time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time) - 1
      
      # Specify years to read full file if no value is provided
      if(is.na(startdate)){startdate <- time[1]}
      if(is.na(enddate)){enddate <- time[length(time)]}
      
      # Set start and count to read specific x and y steps.
      start[1] <- min(which(x >= terra::xmin(ext) & x <= terra::xmax(ext)))
      start[2] <- min(which(y >= terra::ymin(ext) & y <= terra::ymax(ext)))
      start <- as.numeric(start)
      count[1] <- varsize[1] - start[1] - (varsize[1] - max(which(x >= terra::xmin(ext) & x <= terra::xmax(ext))))
      count[2] <- varsize[2] - start[2] - (varsize[2] - max(which(y >= terra::ymin(ext) & y <= terra::ymax(ext))))
      
      # Specify which dates to get
      if(any(time >= startdate) & any(time <= enddate)){
        # Select start and count value to correspond with start and enddate
        start[ndims] <- min(which(time >= startdate & time <= enddate))
        start <- as.numeric(start)
        end <- max(which(time >= startdate & time <= enddate))
        count[ndims] <- end - start[ndims] + 1
        time <- time[start[ndims]:end]
        
        # Read all values for subset of coordinates
        z <- ncdf4::ncvar_get(nc, start=start, count=count)
        ncdf4::nc_close(nc)
        
        # Convert array into raster stack
        z <- terra::rast(z, extent=c(xmin(ext), xmax(ext), ymin(ext),ymax(ext)), 
                         crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
        names(z) <- time
        terra::time(z) <- time
        return(z)
      }
    })
    
    # Remove lists with Null value
    data <- Filter(Negate(is.null), data)
    
    # Turn into one stack
    data <- terra::rast(data)
    
    if(is(mask, "SpatVector")){
      data <- terra::mask(data, mask)
    }
    # Save to file
    if(filename != ""){
      terra::writeCDF(data, filename=filename, overwrite=overwrite, zname="Date")
    }
  }
  return(data)
}