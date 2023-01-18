#' Summarise values of one or multiple NetCDF Files over time and/or space
#'
#' Summarise spatial/spatio-temporal measurements of one or multiple NetCDF files 
#' to a raster stack or data.frame with desired spatial and temporal resolution
#' 
#' @param files \code{character}. A filepath or list of filepaths. Filepath should lead to a NetCDF file.
#' @param startdate \code{integer}. Start year.
#' @param enddate \code{integer}. End year.
#' @param stat \code{character}. Currently only works with mean.
#' @param ext \code{extent}. If present the NetCDF file is subset by this extent.
#' @param cores \code{integer}. Number of cores for parallel computing. If this argument is not provided, at least half the number of availabe cores will be used.
#' @param filename \code{character}. Output filename of the averaged data. If this argument is not provided, result will not be written to disk.
#' @return A dataframe with temporal values averaged over the specified area.
#' @examples
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"),
#'                     pattern="tas.*\\.nc", full.names=TRUE)
#' cellstatsNC(files, startdate=2000, enddate=NA, cores=1)
#' @export cellstatsNC
#' @name cellstatsNC
cellstatsNC <- function(files, startdate=NA, enddate=NA, ext="", cores=NA, filename='', stat="mean"){
  if(filename!="" & file.exists(filename)){
    r_z <- utils::read.csv(filename)
  } else{
    mask <- NA
    if(is(ext,"Extent")){
      # Extent is already an extent object
    } else if(is(ext, "SpatialPolygonsDataFrame")){
      mask <- ext
      ext <- raster::extent(ext)
    } else if(is(ext, "RasterLayer")){
      mask <- ext
      ext <- raster::extent(ext)
    } else if(any(ext != "")){
      ext <- raster::extent(ext)
    }
    
    # Obtain size and other specifications of file
    nc <- ncdf4::nc_open(files[1])
    var <- nc$var[[nc$nvars]]
    varsize <- var$varsize
    ndims   <- var$ndims
    
    # Obtain coordinates
    lon <- ncdf4::ncvar_get(nc, var$dim[[1]]$name, start=1)
    lat <- ncdf4::ncvar_get(nc, var$dim[[2]]$name, start=1)
    
    #Set count and start to	1,1,1,...,1
    start <- rep(1,ndims)
    count <- rep(1,ndims)
    
    # Specify x, y coordinates
    if(is(ext, "Extent")){
      start[1] <- min(which(lon >= raster::xmin(ext) & lon <=  raster::xmax(ext)))
      start <- as.numeric(start)
      end <-  max(which(lon >=  raster::xmin(ext) & lon <=  raster::xmax(ext)))
      count[1] <- end - start[1] + 1
      start[2] <- min(which(lat >=  raster::ymin(ext) & lat <=  raster::ymax(ext)))
      start <- as.numeric(start)
      end <-  max(which(lat >=  raster::ymin(ext) & lat <=  raster::ymax(ext)))
      count[2] <- end - start[1] + 1
    } else{
      count[1] <- varsize[1]
      count[2] <- varsize[2]
    }
    
    #Get variable name of file
    name <- var$name
    
    timeref <- as.Date(strsplit(nc$dim[[nc$ndims]]$units, " ")[[1]][3]) 
    if(ncdf4::ncvar_get(nc, nc$dim$time)[1] == 0){
      time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
    } else if(ncdf4::ncvar_get(nc, nc$dim$time)[1] == 1){
      time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time) - 1
    }
    
    # Define start date
    if(is.na(startdate)){
      startdate <- time[1]
    } else if(is(startdate, "Date")){
      # startdate is already in Date format
    } else {
      startdate <- as.Date(paste0(startdate, "-01-01"))
    }
    
    # Close NC file again
    ncdf4::nc_close(nc)
    
    # Define end date
    if(is.na(enddate)){
      nc <- ncdf4::nc_open(files[length(files)])
      timeref <- as.Date(strsplit(nc$dim[[nc$ndims]]$units, " ")[[1]][3]) 
      if(ncdf4::ncvar_get(nc, nc$dim$time)[1] == 0){
        time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
      }else if (ncdf4::ncvar_get(nc, nc$dim$time)[1] == 1){
        time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time) - 1
      }
      enddate <- time[length(time)]
      ncdf4::nc_close(nc)
    } else if(is(enddate, "Date")){
      # enddate is already in date formate
    } else{
      enddate <-  as.Date(paste0(enddate, "-12-31"))
    }
    
    # Calculate the number of cores available and use 0.55
    if(is.na(cores)){cores <- ceiling(0.55*parallel::detectCores())}
    
    # Initiate cluster
    cl <- parallel::makeCluster(cores)
    
    # Load variables
    parallel::clusterExport(cl, list("files", "start", "count", "name", "cores", 
                                     "stat", "startdate", "enddate"), 
                            envir=environment())
    
    # Load packages for cluster
    parallel::clusterEvalQ(cl, library(ncdf4))
    parallel::clusterEvalQ(cl, library(abind))
    parallel::clusterEvalQ(cl, library(raster))
    
    #Run in parallel
    data <- parallel::parLapply(cl, files, function(a){
      data <- lapply(1:4, FUN=function(z){
        nc <- ncdf4::nc_open(a)
        var <- nc$var[[nc$nvars]]
        varsize <- var$varsize
        ndims   <- var$ndims
        
        # Get start and end date of file
        timeval <- ncdf4::ncvar_get(nc, var$dim[[ndims]]$name, start=1, 
                                    count=varsize[ndims])
        timeref <- as.Date(strsplit(var$dim[[ndims]]$units, " ")[[1]][3]) 
        time <- timeref + timeval
        
        if(any(time >= startdate) & any(time <= enddate)){
          # Select start and count value to correspond with start and enddate
          start[ndims] <- min(which(time >= startdate & time <= enddate))
          start <- as.numeric(start)
          end <- max(which(time >= startdate & time <= enddate))
          count[ndims] <- end - start[ndims] + 1
          time <- time[start[ndims]:end]
          
          # Split data by 4
          splits <- round(seq(0, count[ndims], length.out=4+1))
          start[ndims] <- start[ndims] + splits[z]
          count[ndims] <- splits[z+1] - splits[z] - 1
          
          # Get data
          data <- ncdf4::ncvar_get(nc, var, start=start, count=count)
          ncdf4::nc_close(nc)
          
          # Mask data and calculate mean
          if(is(mask, "SpatialPolygonsDataFrame") | is(mask, "RasterLayer")){
            # Re-arrange array
            data <- aperm(data, c(2,1,3))
            
            # Turn matrix into raster brick
            data <- raster::brick(data, xmn=raster::xmin(ext), xmx=raster::xmax(ext), 
                                  ymn=raster::ymin(ext), ymx=raster::ymax(ext), 
                                  crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
            # Mask data
            data <- raster::mask(data, mask)
            
            # Calculate mean
            data <- as.data.frame(raster::cellStats(data, stat=stat)); gc()
          } else if(stat == "mean"){
            # Calculate mean
            data <- as.data.frame(apply(data, 3, FUN=function(x) mean(x, na.rm = TRUE))); gc()
          }
          return(data)
        }
      })
      # Remove lists with Null value
      data <- Filter(Negate(is.null), data)
      
      # Combine matrices
      data <- do.call("rbind", data)
      
      return(data)
    })
    
    # Close the cluster
    parallel::stopCluster(cl)
    
    # Combine matrices
    data <- do.call("rbind", data)
    colnames(data) <- stat
    data$date <- seq(startdate, enddate, length.out=nrow(data))
    
    # Save to file
    if(filename != ""){
      utils::write.csv(data, path=filename)
    }
  }
  return(data)
}