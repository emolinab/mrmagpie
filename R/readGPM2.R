#' @title readGPM2
#' @description read peatland area from GPM2
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#' \dontrun{
#' readSource("GPM2", convert = "onlycorrect")
#' }
#' @importFrom magclass as.magpie

readGPM2 <- function() {
  terra::terraOptions(tempdir = withr::local_tempdir(tmpdir = getConfig("tmpfolder")),
                      todisk = FALSE, memfrac = 0.5)
  withr::defer(terra::terraOptions(tempdir = tempdir()))

  # read-in file
  r <- terra::rast("peatMAY22_1x1_mw_RUS30.tif")

  # choose only 1st layer
  r1 <- terra::segregate(r, other = NA)
  r1 <- r1[[1]]

  # project r to lon/lat
  r2 <- terra::project(r1, "+proj=longlat +datum=WGS84", method = "near")

  # get cell area
  a <- terra::cellSize(r2, unit = "ha", mask = TRUE)
  a <- a * 1e-6

  # project or aggregate to 0.5 degree
  # use terra::aggregate because terra::project(a, terra::rast(res = 0.5), method = "sum") is not working
  # on the cluster (method = "sum" is the problem)
  r3 <- terra::aggregate(a, fact = 48, fun = "sum", na.rm = TRUE)

  # get spatial mapping
  map <- mrcommons::toolGetMappingCoord2Country(pretty = TRUE)

  # transform raster to magpie object
  x <- as.magpie(terra::extract(r3, map[c("lon", "lat")])[, -1], spatial = 1)

  # set dimension names
  dimnames(x) <- list("coords" = map$coords, "t" = NULL, "d3" = NULL)

  return(x)
}
