#' @title calcExtremesFrequency
#' @description It calculates the number of extreams events occurring in a define number of years.
#' E.g., Occurence over 20 years.
#' @param subtype details of climate model, scenario, and extreme variable in the format: gcm:scenario:variable
#' @param initialYear Initial year of the magpie object, e.g., 1995
#' @param yearsOver Occurence over xx years. E.g., Occurence over 20 years
#' @param gcmAvg when all climate models are read, should the funcion average overall of them
#' @param percentage should the output be given as a percentage (TRUE) or years over xx years (FALSE) or
#' as fraction ("fraction")
#' @param method the source of the extremes data
#' @return magpie object
#' @author Edna J. Molina Bacca
#' @examples
#' \dontrun{
#' calcOutput("ExtremesFrequency", subtype = "MRI-ESM2-0:ssp370:multi", initialYear = 1995)
#' }
#'
#' @importFrom magclass as.magpie

calcExtremesFrequency <- function(subtype = "MRI-ESM2-0:ssp370:multi", initialYear = 1995,
                                  yearsOver = 30, gcmAvg = FALSE, percentage = FALSE, method = "Biess2024") {


  if (method == "Biess2024") {
    x <- readSource(method, subtype = subtype, subset = seq(initialYear - (yearsOver - 1), 2100, 1))
  } else {
    stop("The source is not available")
  }

  subtype <- toolSplitSubtype(subtype,
                              list(climatemodel = NULL,
                                   scenario     = NULL,
                                   variable     = NULL))

  yrs <- as.numeric(sub("y", "", getYears(x)))
  targets <- yrs[yrs >= initialYear]

  if (subtype$climatemodel == "all") {
    if (isTRUE(gcmAvg)) {
      x <- round(magpply(x, mean, MARGIN = c(1, 2, 3.2, 3.3), na.rm = TRUE), 0)
      getNames(x) <- paste0("avg.", getNames(x))
    } else if (!isFALSE(gcmAvg)) {
      stop("The gcmAvg argument must be TRUE or FALSE")
    }
  }

  neededYears <- paste0("y", seq(min(targets) - (yearsOver - 1), max(targets)))
  x <- x[, neededYears, , drop = FALSE]

  outList <- lapply(targets, function(y) {
    xSub <- x[, paste0("y", (y - (yearsOver - 1)):y), , drop = FALSE]
    xSum <- magpply(xSub, sum, MARGIN = c(1, 3), na.rm = TRUE)
    getYears(xSum) <- paste0("y", y)
    xSum
  })
  out <- mbind(outList)

  out <- switch(as.character(percentage),
    "fraction" = out / yearsOver,
    "TRUE"     = out / yearsOver * 100,
    "FALSE"    = out
  )

  unit <- switch(as.character(percentage),
    "fraction" = "Fraction",
    "TRUE"     = "Percentage",
    "FALSE"    = paste0("Number of years out of ", yearsOver)
  )
  weight <-  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,
                        aggregate = FALSE)[, , "q33_marginal"]

  return(list(
    x = out,
    weight = weight,
    unit = unit,
    description = "Frequency in the occurence of extreme events",
    isocountries = FALSE
  ))
}
