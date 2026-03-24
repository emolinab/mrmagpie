#' @title calcYieldAnomaly
#' @description It calculates the anomalies from a [certain number of years] moving average for LPJmL data,
#' calculated as (x-xAvg)/xAvg
#' @param subtype details of climate model and scenario
#' @param initialYear Initial year of the magpie object, e.g., 1995
#' @param endYear End year of the magpie object, e.g., 2100
#' @param yearsOver Window of years for the moving average. E.g., 30
#' @param kcr "all" or specific crop type to calculate vulnerability
#' @param version "version of LPJmL to use
#' @param timeStep years difference between time steps of the output
#' @return magpie object
#' @author Edna J. Molina Bacca
#' @examples
#' \dontrun{
#' calcOutput("calcYieldAnomaly", subtype = "MRI-ESM2-0:ssp370", initialYear = 1995,
#'                                  yearsOver = 20, kcr = "all")
#' }
#'
#' @importFrom magclass as.magpie
#' @importFrom madrat toolSplitSubtype
#' @importFrom data.table as.data.table  frollmean
#' @importFrom data.table :=

calcYieldAnomaly <- function(subtype = "MRI-ESM2-0:ssp370", initialYear = 1995, endYear = 2100,
                             yearsOver = 30, kcr = "all", version = "ggcmi_phase3_nchecks_bft_6277d36e",
                             timeStep = 5) {

  year <- celliso <- data <- value <- rollingMean <- anomaly <- NULL

  # LPJmL raw data from the initialYear - yearsOver + 1
  lpjmLData <- calcOutput("LPJmL_new", version = version, climatetype = subtype,
                          subtype = "harvest", stage = "raw", aggregate = FALSE)[, seq((initialYear - yearsOver + 1),
                                                                                       endYear, 1), ]

  # data table to be able to manage the large set
  lpjmLData <- as.data.table(lpjmLData, spatial = TRUE, temporal = TRUE)
  setnames(lpjmLData, c("x.y.iso", "data.data1", "Year"), c("celliso", "data", "year"))
  lpjmLData <- lpjmLData[, year := as.integer(sub("y", "", year))][, list(celliso, data, year, value)] # nolint: object_usage_linter

  lpjmLData <-  lpjmLData[, rollingMean := frollmean(value, n = yearsOver, align = "right"), # nolint: object_usage_linter
                        by = list(celliso, data)][, list(celliso, data, year, value, rollingMean)] # nolint: object_usage_linter

  # Anomaly as the relative difference between value and average
  lpjmLData <-  lpjmLData[, `:=`(value = round(value, 2), # nolint: object_usage_linter
                                 anomaly = round((value - rollingMean) / rollingMean, 1))][, # nolint: object_usage_linter
                               list(celliso, data, year, anomaly)][year %in% seq(initialYear, endYear, timeStep)] # nolint: object_usage_linter

  out <- as.magpie(lpjmLData, spatial = "celliso", temporal = "year")
  getCells(out) <- gsub("_", "\\.", getCells(out))
  getNames(out) <- gsub("_", "\\.", getNames(out))

  weight <-  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
                        aggregate = FALSE)[, , "q33_marginal"]

  return(list(
    x = out,
    weight = weight,
    unit = paste0("Fraction"),
    description = "Relative difference from a", yearsOver, " moving average",
    isocountries = FALSE
  ))
}
