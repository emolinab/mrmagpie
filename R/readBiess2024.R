#' @title readBiess2024
#' @description Biess data contains information on 2.5° x 2.5° on occurence of extreme events based on different
#'              CMIP6 gcm data. This information based on the paper doi: 10.1088/1748-9326/ad4619.
#' @param subtype climatemodel:scenario:variable
#'           - climate model: specific (e.g., "MRI-ESM2-0") or "all"
#'           - scenario: ssp126, ssp245, ssp370 available
#'           - variable: Occurrence of "p_e" (dry) , "pr" (wet), or "tasmax" (heat) extreams or multi
#'             (any of the three) (0 or 1) in one year
#'
#' @param subset Years to read as a vector
#' @return MAgPIE object of the defined data
#' @author Edna Molina Bacca
#' @examples
#' \dontrun{
#' readSource("Biess2024", subtype = "MRI-ESM2-0:ssp370:multi", subset = seq(1995, 2100, 1))
#' }
#' @importFrom madrat toolSplitSubtype
#' @importFrom terra rast varnames disagg
#' @importFrom mstools toolGetMappingCoord2Country
#' @importFrom dplyr %>% mutate left_join filter select if_else
#' @importFrom tidyr pivot_longer

readBiess2024 <- function(subtype = "MRI-ESM2-0:ssp370:multi", subset = seq(1966, 2100, 1)) {


  subtype <- toolSplitSubtype(subtype,
                              list(climatemodel = NULL,
                                   scenario     = NULL,
                                   variable     = NULL))

  mapping <- toolGetMappingCoord2Country(pretty = TRUE, extended = FALSE)

  # Reads the specific file
  rasterMulti <- rast(paste0(subtype$variable, "_", subtype$scenario, ".nc"))

  # Variable names
  varNames <- varnames(rasterMulti)
  varNames <- if (subtype$climatemodel == "all") varNames else grep(subtype$climatemodel, varNames, value = TRUE)

  # Subset of the data
  years <- subset - 1849
  namT <- paste0(varNames, "_", years)
  rasterMulti <- rasterMulti[[namT]]

  rasterMulti <- disagg(rasterMulti, fact = 5)

  rasterMAg <- as.data.frame(rasterMulti, xy = TRUE)


  suppressWarnings({
    rasterMAg1 <- rasterMAg |>
      pivot_longer(
        cols = -c("x", "y"),                 # All columns but lat and lon
        names_to = c("model", "year"),   # Split model and year
        names_pattern = "(.*)_(\\d+)",   # Split by second underscore
        values_to = "Value"
      ) |>
      mutate(year = as.integer(.data$year) + 1849,
             lon   = if_else(x > 180, x - 360, x),  # Corrects longitudes >180
             lat   = .data$y,
             Data2 = subtype$scenario,
             Data3 = subtype$variable,
             Data1 = paste(sub("_.*", "", .data$model), .data$Data2, .data$Data3, sep = "."), )  |>
      left_join(mapping |>
                  select(.data$iso, .data$lon, .data$lat), by = c("lon", "lat")) |>
      mutate(
        lon      = if_else(.data$x > 180, .data$x - 360, .data$x),
        lon      = gsub("\\.", "p", .data$lon),
        lat      = gsub("\\.", "p", .data$lat),
        celliso  = paste0(.data$lon, ".", .data$lat, ".", .data$iso)
      ) |>
      filter(!is.na(.data$iso)) |>
      select(.data$celliso, Year = .data$year, .data$Data1, .data$Value) |>
      as.data.frame()
  })


  x <- as.magpie(rasterMAg1, temporal = "Year", spatial = "celliso")
  getNames(x) <- gsub("_", "\\.", getNames(x))
  getCells(x) <- gsub("_", "\\.", getCells(x))
  x <- magpiesort(x)


  return(x)
}
