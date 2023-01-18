#' R Package for processing large NetCDF files
#'
#' Analysis of large NetCDF files (Cropping, Summarising) for 
#' further use with R using the raster package.
#' 
#' @name processNC-package
#' @aliases processNCpackage
#' @docType package
#' @title R Package for processing large NetCDF files
#' @author RS-eco
#'
#' @import ncdf4 raster sp parallel rlang
#' @importFrom tidyselect one_of
#' @importFrom magrittr %>%
#' @importFrom abind abind
#' @importFrom tidyr gather
#' @importFrom dplyr funs group_by_at group_by summarise_all vars
#' @importFrom lubridate year month
#' @importFrom utils globalVariables
#' @importFrom methods is
#' @keywords package
#'
NULL
#'
#' @docType data
#' @name bavaria
#' @title Global Administrative Boundary of Bavaria
#' @description Global Administrative Boundary (GADM) of Bavaria
#' @details This dataset is a shapefile containing a single polygon 
#' of the administrative boundary of Bavaria.
#' @format \code{sp::SpatialPolygonsDataFrame}
#' 
NULL