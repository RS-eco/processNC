#' ## Create example NetCDF file

# Data can be downloaded from: https://data.isimip.org/10.5880/pik.2019.004

var <- "pr" # pr, tas

# Get daily temperature EWEMBI ISIMIP2b files
files <- list.files("/home/matt/Documents/EWEMBI", 
                        pattern=paste0("^", var, "_"), full.names=T)

# Load data
data <- lapply(files, raster::stack)

# Crop data by extent of Germany
deu <- raster::getData("GADM", country="DEU", level=1, path="/home/matt/Documents/")
dat_deu <- lapply(data, function(x) raster::mask(raster::crop(x, deu), deu))

# Get list of z units
z_units <- lapply(files, function(x) ncdf4::nc_open(x)$dim[[3]]$units)

# Save data to NetCDF file
mapply(FUN=function(x,y,z){
  if(var == "tas"){
    raster::writeRaster(x=x, filename=paste0("inst/extdata/", y), format="CDF",
           varname="tas", varunit="kg m-2 s-1", xname="lon", yname="lat", 
           zname="time", zunit=z, overwrite=T)
  } else if(var == "pr"){
    raster::writeRaster(x=x, filename=paste0("inst/extdata/", y), format="CDF",
                        varname="pr", varunit="Kelvin", xname="lon", yname="lat", 
                        zname="time", zunit=z, overwrite=T)
  }
}, x=dat_deu, y=sub("ewembi", "ewembi_deu", basename(files)), z=z_units)

files <- paste0("inst/extdata/", sub(".nc4", ".nc", sub("ewembi", "ewembi_deu", 
                                                         basename(files))))
nc <- ncdf4::nc_open(files[1])
nc$dim[[3]]$units
