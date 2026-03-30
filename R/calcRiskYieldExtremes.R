#' @title calcRiskYieldExtremes
#' @description It calculates the anomalies from a [certain number of years] moving average for LPJmL data,
#' calculated as (x-xAvg)/xAvg
#' @param subtype details of climate model and scenario
#' @param initialYear Initial year of the magpie object. Due to data availability first year available is 2010
#' @param endYear End year of the magpie object, e.g., 2100
#' @param yearsOver Window of years for the moving average. E.g., 30
#' @param kcr "all" or specific crop type to calculate vulnerability
#' @param version "version of LPJmL to use
#' @param timeStep years difference between time steps of the output
#' @param extremeType Type of extreme "p_e" (dry) , "pr" (wet), or "tasmax" (heat) extreams or multi
#'             (any of the three)
#' @return magpie object
#' @author Edna J. Molina Bacca
#' @examples
#' \dontrun{
#' calcOutput("calcRiskYieldExtremes", subtype = "MRI-ESM2-0:ssp370", initialYear = 1995,
#'                                  yearsOver = 20, kcr = "all")
#' }
#'
#' @importFrom magclass as.magpie
#' @importFrom madrat toolSplitSubtype
#' @import data.table
#' @importFrom stats wilcox.test

calcRiskYieldExtremes <- function(subtype = "MRI-ESM2-0:ssp370", initialYear = 2025, endYear = 2100,
                                  yearsOver = 30, kcr = "all", version = "ggcmi_phase3_nchecks_bft_6277d36e",
                                  timeStep = 5, extremeType = "multi") {

  celliso <- year <- value <- year <- valueHazard <- data <- celliso <- eHaz <- eNoHaz <- NULL

  # As data table to speed the calculation up
  yieldAnomaly <- as.data.table(calcOutput("YieldAnomaly", aggregate = FALSE, subtype = subtype,
                                           initialYear = (initialYear - yearsOver + 1), endYear = endYear,
                                           yearsOver = yearsOver, kcr = kcr, version = version, timeStep = 1))

  setnames(yieldAnomaly, c("x.y.iso", "data.data1"), c("celliso", "data"))

  binaryHazard <- as.data.table(readSource("Biess2024", subtype = paste0(subtype, ":", extremeType),
                                           subset = seq((initialYear - yearsOver + 1), endYear, 1)))
  setnames(binaryHazard, c("celliso.region.region1", "Year", "Data1.data.data1"), c("celliso", "year", "data"))

  anomHazDT <- merge(
    yieldAnomaly,
    binaryHazard[, list(celliso, year, value)],
    by = c("celliso", "year"),
    suffixes = c("", "Hazard")
  )

  anomHazDT <- anomHazDT[, `:=`(year = as.numeric(gsub("y", "", year)))]

  expectFunction <- function(x) {
    freq <- table(round(x, 0))
    sum(as.numeric(names(freq)) * freq) / sum(freq)
  }
  targetYears <- seq(initialYear, endYear, by = timeStep)


  resultList <- lapply(targetYears, function(y) {
    dfWindow <- anomHazDT[year <= y & year > (y - yearsOver + 1)]
    dfWindow[, list(
      eHaz = expectFunction(value[valueHazard == 1]),
      eNoHaz = expectFunction(value[valueHazard == 0])
    ), by = list(celliso, data)][, year := y]
  })

  expectHaz <- rbindlist(resultList)
  expectHaz <- expectHaz[eNoHaz < eHaz, eHaz := 0]

  # Risk from the mean

  hazProbability <- as.data.table(calcOutput("ExtremesFrequency", subtype = paste0(subtype, ":", extremeType),
                                             initialYear = initialYear,
                                             yearsOver = yearsOver, percentage = "fraction", method = "Biess2024",
                                             aggregate = FALSE))
  setnames(hazProbability, c("x.y.iso", "Year"), c("celliso", "year"))
  hazProbability <- hazProbability[, year := as.integer(sub("y", "", year))][, list(celliso, year, value)]

  riskDT <- merge(hazProbability[year %in% targetYears],
                  expectHaz[year %in% targetYears][, list(celliso, year, data, eHaz, value)],
                  by = c("celliso", "year"))[, `:=`(value = value * eHaz)
  ][, list(celliso, year, data,
           value)][value > 0 | !(is.finite(value)), `:=`(value = 0)]

  out <- as.magpie(riskDT, spatial = "celliso", temporal = "year")
  getNames(out) <- gsub("_", "\\.", getNames(out))
  getCells(out) <- gsub("_", "\\.", getCells(out))
  out[!is.finite(out)] <- 0


  # Make the risk zero for those years-crops combination were pValue is larger than 0.05
  # (meaning that the Haz - noHaz data is statistically are not different)

  pValue <- anomHazDT[year %in% targetYears][, list(
                                                    pValue = wilcox.test(
                                                      .SD[valueHazard == 1]$value, # & value!=0
                                                      .SD[valueHazard == 0]$value
                                                    )$p.value), by = list(data, year)]


  nonSig <- pValue[pValue > 0.05 | is.na(pValue)]

  for (i in seq_len(nrow(nonSig))) {
    out[, nonSig$year[i], nonSig$data[i]] <- 0
  }


  # Mapping LPJmL to MAgPIE crops
  lpj2mag   <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")
  out <- toolAggregate(out, lpj2mag, from = "LPJmL",
                       to = "MAgPIE", dim = 3.1, partrel = TRUE)
  out <- out[, , "pasture", invert = TRUE]

  weight <-  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
                        aggregate = FALSE)[, , "q33_marginal"]


  return(list(
    x = out,
    weight = weight,
    unit = paste0("Fraction of yield at risk"),
    description = "Fraction of yield at risk due to extreme events based on a probabilistic risk analysis",
    isocountries = FALSE
  ))
}
