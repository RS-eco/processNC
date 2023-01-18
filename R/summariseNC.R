#' Summarise values of one or multiple NetCDF Files over time and/or space
#'
#' Summarise spatial/spatio-temporal measurements of one or multiple NetCDF files 
#' to a raster stack or data.frame with desired spatial and temporal resolution
#' 
#' @param files \code{character}. A filepath or list of filepaths. Filepath should lead to a NetCDF file.
#' @param startdate \code{integer}. Start year.
#' @param enddate \code{integer}. End year.
#' @param group_col \code{character}. Columns to use for group_at call, one of week, month, year, c("week", "year") or c("month", "year")
#' @param ext \code{extent}. If present the NetCDF file is subset by this extent.
#' @param name \code{character}. Name of variable, can be one of pr, tas, tasmax, tasmin, sfcWind, hurs, huss, dis
#' @param cores \code{integer}. Number of cores for parallel computing. If this argument is not provided, at least half the number of availabe cores will be used.
#' @param filename1 \code{character}. Output filename of the averaged data. If this argument is not provided, result will not be written to disk.
#' @param filename2 \code{character}. Output filename of the coefficient of variation. If this argument is not provided, the coefficient of variation will not be written to disk.
#' @param overwrite \code{logical}. If TRUE, existing files will be overwritten.
#' @return A raster stack with monthly layers of aggregated data over the specified time period and area.
#' @examples
#' files <- list.files(paste0(system.file(package="processNC"), "/extdata"), 
#'                     pattern="tas.*\\.nc", full.names=TRUE)
#' summariseNC(files=files[4], group_col=c("year", "month"), 
#'             startdate=2001, enddate=2010, cores=1)
#' summariseNC(files=files[4], group_col=c("year", "week"), 
#'             startdate=2001, enddate=2010, cores=1)
#' summariseNC(files=files[4], group_col="year", 
#'             startdate=2001, enddate=2010, cores=1)
#' summariseNC(files=files[4], group_col="month", 
#'             startdate=2001, enddate=2010, cores=1)
#' summariseNC(files=files[4], group_col="week", 
#'             startdate=2001, enddate=2010, cores=1)
#' @seealso [aggregateNC]
#' @export summariseNC
#' @name summariseNC
summariseNC <- function(files, startdate=NA, enddate=NA, ext=NA, group_col=c("year", "month"),
                        name=NA, cores=NA, filename1='', filename2='', overwrite=FALSE){
  . <- NULL; x <- NULL; week <- NULL
  if(overwrite==FALSE & filename1!="" & file.exists(filename1)){
    r_z <- terra::rast(filename1)
  } else{
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
    
    # Obtain size and other specifications of file
    nc <- ncdf4::nc_open(files[1])
    var <- nc$var[[nc$nvars]]
    varsize <- var$varsize
    ndims   <- var$ndims
    
    # Obtain coordinates
    lon <- ncdf4::ncvar_get(nc, var$dim[[1]]$name, start=1)
    lat <- ncdf4::ncvar_get(nc, var$dim[[2]]$name, start=1)
    
    # Specify range of y values
    if(is(ext, "Extent")){
      ny <- which(lat >= terra::ymin(ext) & lat <= terra::ymax(ext))
      lat <- lat[ny]
    } else{
      ny <- seq(1, varsize[2], 1)
      lat <- lat[ny]
    }  
    
    #Get variable name of file
    if(is.na(name)){
      name <- var$name
    }
    
    #Set count and start to	1,1,1,...,1
    count <- rep(1,ndims)
    start <- rep(1,ndims)
    
    # Specify longitude dimensions (either all or specific extent)
    count[1] <- varsize[1]
    if(is(ext, "Extent")){
      start[1] <- min(which(lon >= terra::xmin(ext) & lon <= terra::xmax(ext)))
      start <- as.numeric(start)
      end <-  max(which(lon >= terra::xmin(ext) & lon <= terra::xmax(ext)))
      count[1] <- end - start[1] + 1
      lon <- lon[start[1]:end]
    }
    
    timeref <- as.Date(strsplit(var$dim[[nc$ndims]]$units, " ")[[1]][3]) 
    if(ncdf4::ncvar_get(nc, nc$dim$time)[1] == 0){
      time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
    } else if (ncdf4::ncvar_get(nc, nc$dim$time)[1] == 1){
      time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time) - 1
    } else if (ncdf4::ncvar_get(nc, nc$dim$time)[1] > 1000){
      time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
    }
    
    # Close NC file again
    ncdf4::nc_close(nc)
    
    # Define start date
    if(is.na(startdate)){
      startdate <- time[1]
    } else if(is(startdate, "Date")){
      # startdate is already in Date format
    } else {
      startdate <- as.Date(paste0(startdate, "-01-01"))
    }
    
    # Define end date
    if(is.na(enddate)){
      nc <- ncdf4::nc_open(files[length(files)])
      timeref <- as.Date(strsplit(nc$dim[[nc$ndims]]$units, " ")[[1]][3]) 
      if(ncdf4::ncvar_get(nc, nc$dim$time)[1] == 0){
        time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
      } else if (ncdf4::ncvar_get(nc, nc$dim$time)[1] != 1){
        time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time) - 1
      }
      enddate <- time[length(time)]
      ncdf4::nc_close(nc)
    } else if(is(enddate, "Date")){
      # enddate is already in Date format
    } else{
      enddate <-  as.Date(paste0(enddate, "-12-31"))
    }
    
    # Calculate the number of cores available and leave one for basic use
    if(is.na(cores)){cores <- ceiling(0.55*parallel::detectCores())}
    
    # Initiate cluster
    cl <- parallel::makeCluster(cores)
    
    # Load variables
    parallel::clusterExport(cl, list("ny", "files", "start", "count", "name", 
                                     "lon", "lat", "group_col", "startdate", "enddate"), 
                            envir=environment())
    
    # Load packages for cluster
    parallel::clusterEvalQ(cl, library(ncdf4))
    parallel::clusterEvalQ(cl, library(lubridate))
    parallel::clusterEvalQ(cl, library(terra))
    parallel::clusterEvalQ(cl, library(dplyr))
    
    #Run in parallel
    z <- parallel::parLapply(cl, ny, function(y){
      data <- lapply(files, FUN=function(a){
        nc <- ncdf4::nc_open(a)
        var <- nc$var[[nc$nvars]]
        varsize <- var$varsize
        ndims   <- var$ndims
        
        # Initialize start to read all x and specific y step of the variable.
        start[2] <- y
        
        # Get start and end date of file
        timeref <- as.Date(strsplit(var$dim[[nc$ndims]]$units, " ")[[1]][3]) 
        if(ncdf4::ncvar_get(nc, nc$dim$time)[1] == 0){
          time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
        } else if (ncdf4::ncvar_get(nc, nc$dim$time)[1] == 1){
          time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time) - 1
        } else if (ncdf4::ncvar_get(nc, nc$dim$time)[1] > 1000){
          time <- timeref + ncdf4::ncvar_get(nc, nc$dim$time)
        }
        
        if(any(time >= startdate) & any(time <= enddate)){
          # Select start and count value to correspond with start and enddate
          start[ndims] <- min(which(time >= startdate & time <= enddate))
          start <- as.numeric(start)
          end <-  max(which(time >= startdate & time <= enddate))
          count[ndims] <- end - start[ndims] + 1
          time <- time[start[ndims]:end]
          
          # Get data
          data <- ncdf4::ncvar_get(nc, var, start=start, count=count)
          ncdf4::nc_close(nc)
          return(data)
        }
      })
      
      data <- as.data.frame(t(do.call("cbind", data)))
      colnames(data) <- lon
      data$date <- seq(startdate, enddate, length.out=nrow(data))
      data$y <- lat[y]
      
      # Require dplyr package
      requireNamespace("dplyr")
      
      # Specify function
      if(name %in% c("hurs", "huss", "tas", "sfcWind")){
        .funs = list(~ mean(., na.rm = TRUE))
      } else if(name %in% c("pr", "dis")){
        .funs = list(~ sum(., na.rm = TRUE))
      } else if(name == "tasmax"){
        .funs <- list(~ max(., na.rm = TRUE))
      } else if(name == "tasmin"){
        .funs <- list(~ min(., na.rm = TRUE))
      } else{
        .funs = list(~ mean(., na.rm = TRUE))
        print("Note as the variable name is none of the ones specified, the mean is calculated per default.")
      }
      
      # Summarise according to group_col
      if(identical(group_col, c("month", "year")) | identical(group_col, c("year", "month"))){
        data$year <- lubridate::year(data$date)
        data$month <- lubridate::month(data$date)
        # Remove date column for summary
        data <- subset(data, select=-date)
        z_y <- dplyr::group_by_at(data, vars(c(group_col, "y")))
        z_y <- as.data.frame(dplyr::summarise_all(z_y, .funs=.funs))
        z_cv <- dplyr::group_by_at(data, vars(c(group_col, "y")))
        z_cv <- as.data.frame(dplyr::summarise_all(z_cv, funs(raster::cv(., na.rm = TRUE))))
        z <- list()
        z_y <- group_by(z_y, month)
        z$avg <- as.data.frame(summarise_all(z_y, funs(mean(., na.rm = TRUE))))
        z_cv <- group_by(z_cv, month)
        z$cv <- as.data.frame(summarise_all(z_cv, funs(mean(., na.rm = TRUE))))
      } else if(identical(group_col, c("week", "year")) | identical(group_col, c("year", "week"))){
        data$year <- lubridate::year(data$date)
        data$week <- lubridate::week(data$date)
        # Remove date column for summary
        data <- subset(data, select=-date)
        z_y <- dplyr::group_by_at(data, vars(c(group_col, "y")))
        z_y <- as.data.frame(dplyr::summarise_all(z_y, .funs=.funs))
        z_cv <- dplyr::group_by_at(data, vars(c(group_col, "y")))
        z_cv <- as.data.frame(dplyr::summarise_all(z_cv, funs(raster::cv(., na.rm = TRUE))))
        z <- list()
        z_y <- group_by(z_y, week)
        z$avg <- as.data.frame(summarise_all(z_y, funs(mean(., na.rm = TRUE))))
        z_cv <- group_by(z_cv, week)
        z$cv <- as.data.frame(summarise_all(z_cv, funs(mean(., na.rm = TRUE))))
      } else if(group_col == "week"){
        data$year <- lubridate::year(data$date)
        data$week <- lubridate::week(data$date)
        data$week <- as.character(paste(data$week, data$year, sep=" "))
        # Remove date column for summary
        data <- subset(data, select=-c(year,date))
        z <- list()
        z_y <- dplyr::group_by_at(data, vars(c("week", "y")))
        z$y <- as.data.frame(dplyr::summarise_all(z_y, .funs=.funs))
        z_cv <- dplyr::group_by_at(data, vars(c("week", "y")))
        z$cv <- as.data.frame(dplyr::summarise_all(z_cv, funs(raster::cv(., na.rm = TRUE))))
      } else if(group_col == "month"){
        data$month <- as.character(zoo::as.yearmon(data$date))
        # Remove date column for summary
        data <- subset(data, select=-date)
        z <- list()
        z_y <- dplyr::group_by_at(data, vars(c("month", "y")))
        z$y <- as.data.frame(dplyr::summarise_all(z_y, .funs=.funs))
        z_cv <- dplyr::group_by_at(data, vars(c("month", "y")))
        z$cv <- as.data.frame(dplyr::summarise_all(z_cv, funs(raster::cv(., na.rm = TRUE))))
      } else if(group_col == "year"){
        data$year <- lubridate::year(data$date)
        # Remove date column for summary
        data <- subset(data, select=-date)
        z <- list()
        z_avg <- group_by_at(data, vars(group_col, "y"))
        z$avg <- as.data.frame(summarise_all(z_avg, .funs=.funs))
        z_cv <- group_by_at(data, vars(group_col, "y")) 
        z$cv <- as.data.frame(summarise_all(z_cv, funs(raster::cv(., na.rm = TRUE))))
      }
      return(z)
    })
    # Close the cluster
    parallel::stopCluster(cl)
    
    # Turn into a dataframe
    data <- do.call("c", z)#; rm(z)
    avg <- seq(1,length(data), by=2)
    cv <- seq(2,length(data), by=2)
    data1 <- do.call("rbind", data[avg])
    data2 <- do.call("rbind", data[cv]); rm(data)
    data1 <- tidyr::gather(data1, x, avg, -tidyselect::one_of(c(group_col, "y")))
    data2 <- tidyr::gather(data2, x, cv, -tidyselect::one_of(c(group_col, "y")))
    data <- cbind(data1, data2$cv); rm(data1, data2)
    colnames(data)[ncol(data)] <- "cv"
    data$x <- as.numeric(data$x)
    data$y <- as.numeric(data$y)
    if(identical(group_col, c("month", "year")) | identical(group_col, c("year", "month"))){
      r_z <- dplyr::select(data, -cv) 
      r_z <- tidyr::spread(r_z, month, avg)
      r_cv <-dplyr::select(data, -avg)
      r_cv <- tidyr::spread(r_cv, month, cv)
      r_z <- dplyr::select(r_z, -year)
      r_cv <- dplyr::select(r_cv, -year)
      r_z <- r_z[,c(2,1,3:ncol(r_z))]
      r_z <- terra::rast(r_z, type="xyz")
      names(r_z) <- month.name
      r_cv <- r_cv[,c(2,1,3:ncol(r_cv))]
      r_cv <- terra::rast(r_cv, type="xyz")
      names(r_cv) <- month.name
      if(is(mask, "SpatVector") | is(mask, "SpatRaster")){
        r_z <- terra::mask(r_z, mask)
        r_cv <- terra::mask(r_cv, mask)
      }
    } else if(identical(group_col, c("week", "year")) | identical(group_col, c("year", "week"))){
      r_z <- dplyr::select(data, -cv) 
      r_z <- tidyr::spread(r_z, week, avg)
      r_cv <-dplyr::select(data, -avg)
      r_cv <- tidyr::spread(r_cv, week, cv)
      r_z <- dplyr::select(r_z, -year)
      r_cv <- dplyr::select(r_cv, -year)
      r_z <- r_z[,c(2,1,3:ncol(r_z))]
      r_z <- terra::rast(r_z, type="xyz")
      names(r_z) <- unique(data$week)
      r_cv <- r_cv[,c(2,1,3:ncol(r_cv))]
      r_cv <- terra::rast(r_cv, type="xyz")
      names(r_cv) <- unique(data$week)
      if(is(mask, "SpatVector") | is(mask, "SpatRaster")){
        r_z <- terra::mask(r_z, mask)
        r_cv <- terra::mask(r_cv, mask)
      }
    } else if(identical(group_col, "year")){
      r_z <- dplyr::select(data, -cv)
      r_z <- tidyr::spread(r_z, year, avg)
      r_z <- r_z[,c(2,1,3:ncol(r_z))]
      r_z <- terra::rast(r_z, type="xyz")
      names(r_z) <- unique(data$year)
      r_cv <- dplyr::select(data, -avg)
      r_cv <- tidyr::spread(r_cv, year, cv)
      r_cv <- r_cv[,c(2,1,3:ncol(r_cv))]
      r_cv <- terra::rast(r_cv, type="xyz")
      names(r_cv) <- unique(data$year)
      if(is(mask, "SpatVector") | is(mask, "SpatRaster")){
        r_z <- terra::mask(r_z, mask)
        r_cv <- terra::mask(r_cv, mask)
      }
    } else if(identical(group_col, "month")){
      r_z <- dplyr::select(data, -cv)
      r_z <- tidyr::spread(r_z, month, avg)
      r_z <- r_z[,c(2,1,3:ncol(r_z))]
      r_z <- terra::rast(r_z, type="xyz")
      names(r_z) <- as.character(unique(data$month))
      r_cv <- dplyr::select(data, -avg)
      r_cv <- tidyr::spread(r_cv, month, cv)
      r_cv <- r_cv[,c(2,1,3:ncol(r_cv))]
      r_cv <- terra::rast(r_cv, type="xyz")
      names(r_cv) <- as.character(unique(data$month))
      if(is(mask, "SpatVector") | is(mask, "SpatRaster")){
        r_z <- terra::mask(r_z, mask)
        r_cv <- terra::mask(r_cv, mask)
      }
    } else if(identical(group_col, "week")){
      r_z <- dplyr::select(data, -cv)
      r_z <- tidyr::spread(r_z, week, avg)
      r_z <- r_z[,c(2,1,3:ncol(r_z))]
      r_z <- terra::rast(r_z, type="xyz")
      names(r_z) <- as.character(unique(data$week))
      r_cv <- dplyr::select(data, -avg)
      r_cv <- tidyr::spread(r_cv, week, cv)
      r_cv <- r_cv[,c(2,1,3:ncol(r_cv))]
      r_cv <- terra::rast(r_cv, type="xyz")
      names(r_cv) <- as.character(unique(data$week))
      if(is(mask, "SpatVector") | is(mask, "SpatRaster")){
        r_z <- terra::mask(r_z, mask)
        r_cv <- terra::mask(r_cv, mask)
      }
    }
    if(filename1 != ""){
      terra::writeCDF(r_z, filename=filename1, overwrite=overwrite)
    }
    if(filename2 != ""){
      terra::writeCDF(r_cv, filename=filename2, overwrite=overwrite)
    }
  }
  return(r_z)
}