#' @title calcRiskYieldExtremes
#' @description It calculates the anomalies from a [certain number of years] moving average for LPJmL data, calculated as (x-xAvg)/xAvg
#' @param subtype details of climate model and scenario
#' @param initialYear Initial year of the magpie object, e.g., 1995
#' @param endYear End year of the magpie object, e.g., 2100
#' @param yearsOver Window of years for the moving average. E.g., 30
#' @param kcr "all" or specific crop type to calculate vulnerability
#' @param version "version of LPJmL to use
#' @param timeStep years difference between time steps of the output
#' @param vulType vulnerability of the system "system" (As difference of expected anomaly when there is no hazard and when there is hazard) or from the meam "mean"  
#' @return magpie object
#' @author Edna J. Molina Bacca
#' @examples
#' \dontrun{
#' calcOutput("calcRiskYieldExtremes", subtype = "MRI-ESM2-0:ssp370", initialYear = 1995,
#'                                  yearsOver = 20, kcr="all")
#' }
#'
#' @importFrom magclass as.magpie
#' @importFrom madrat toolSplitSubtype
#' @importFrom data.table as.data.table

calcRiskYieldExtremes <- function(subtype = "MRI-ESM2-0:ssp370", initialYear = 1995, endYear = 2100,
                                  yearsOver = 30, kcr="all",version="ggcmi_phase3_nchecks_bft_6277d36e", timeStep= 5, extremeType= "multi") {

  # As data table to speed the calculation up
  YieldAnomaly <- as.data.table(calcOutput("YieldAnomaly", subtype = subtype, initialYear = (initialYear-yearsOver+1), endYear = endYear,
                                  yearsOver = yearsOver, kcr=kcr,version=version, timeStep= 1))
 
  BinaryHazard <- as.data.table(readSource("Biess2024",subtype = paste0(subtype,":",extremeType), subset = seq((initialYear-yearsOver+1), endYear, 1)))
  setnames(BinaryHazard,"celliso","x.y.iso")

  AnomHazDT <- merge(
  YieldAnomaly,
  BinaryHazard[,.(x.y.iso,Year,value)],#[Year %in% unique(YieldAnomaly$Year)],
  by = c("x.y.iso", "Year"),
  suffixes = c("", "_hazard")
) 
  
  AnomHazDT <- AnomHazDT[, `:=`(Year = as.numeric(gsub("y","",Year)))]

  ExpectFunction <- function(x) {
    freq <- table(round(x, 0))
    sum(as.numeric(names(freq)) * freq) / sum(freq)
  }
  targetYears <- seq(initialYear, endYear, by = 5)
  

  resultList <- lapply(targetYears, function(y) {
  df_window <- AnomHazDT[Year <= y & Year > (y - yearsOver+ 1)]
  df_window[, .(
    EHaz = ExpectFunction(value[value_hazard == 1]),
    ENoHaz = ExpectFunction(value[value_hazard == 0])
  ), by = .(x.y.iso, data.data1)][, Year := y]
  })

   ExpectHaz <- rbindlist(resultList)
 

 # Risk from the mean

   HazProbability <- as.data.table(calcOutput("ExtremesFrequency",subtype = subtype, initialYear = initialYear,
                                  yearsOver = yearsOver, percentage = "fraction", source = "Biess2024"))
    
   HazProbability <- HazProbability[, d2 := as.integer(sub("y","",d2))][,.(celliso,d2,value)]
   setnames(HazProbability ,c("celliso","d2"),c("x.y.iso","Year"))

  RiskDT <- merge(
  HazProbability ,
  ExpectHaz,
  by = c("x.y.iso", "Year"))[,`:=` (Value=value*EHaz)][,.(x.y.iso, Year,data.data1,Value)][Value > 0 | !(is.finite(Value)), `:=`(Value=0)]
  
  out<-as.magpie(RiskDT, spatial= "x.y.iso", temporal="Year")
  getNames(out) <- gsub("_", "\\.", getNames(out))
  getCells(out) <- gsub("_", "\\.", getCells(out))
  out[!is.finite(out)]<-0

  # Mapping LPJmL to MAgPIE crops
  lpj2mag   <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")
  out <- toolAggregate(out, lpj2mag, from = "LPJmL",
                               to = "MAgPIE", dim = 3.1, partrel = TRUE)
  out <- out[,,"pasture",invert=TRUE]                             

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