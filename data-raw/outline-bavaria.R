#' ## Create SpatialPolygonsDataFrame of Bavaria

#' Get GADM data of Germany
deu <- raster::getData("GADM", country="DEU", level=1, path="E:/Data/GADM")

#' Subset data by NAME_1
bavaria <- deu[deu$NAME_1 == "Bayern",]
plot(bavaria)

#' Save to file
save(bavaria, file="data/bavaria.rda", compress="xz")
