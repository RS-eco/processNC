# Replace data file with this example code:
  
# 365 values of temperature, 1 for each day in the year
tas=c(2.766,3.491,4.494,5.912,6.989,7.742,7.919,7.027,5.369,3.562,2.814,2.179)
# RasterBrick with 365 layers based on tas + noise
r <- raster(nc=10, nr=10)
s <- brick(lapply(1:12, function(x) setValues(r, G0dm[x]+100*rnorm(ncell(r)))))
  
# Make a netCDF file, also see saving of data files (see summedProbabilities.R)
  
# Make a simple netCDF file
filename <- "atttest_types.nc"
dim <- ncdim_def( "X", "inches", 1:12 )
var <- ncvar_def( "Data", "unitless", dim, -1 ) 
ncnew <- nc_create(filename, var)
  