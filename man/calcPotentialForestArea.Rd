% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPotentialForestArea.R
\name{calcPotentialForestArea}
\alias{calcPotentialForestArea}
\title{calcPotentialForestArea}
\usage{
calcPotentialForestArea(
  refData = "lpj",
  countryLevel = FALSE,
  cells = "lpjcell",
  lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de"),
  climatetype = "MRI-ESM2-0:ssp370"
)
}
\arguments{
\item{refData}{Determines the reference data that the estimated potential
forest area is derived from (currently only "lpj")}

\item{countryLevel}{Whether output shall be at country level.
Requires aggregate=FALSE in calcOutput.}

\item{cells}{magpiecell (59199 cells) or lpjcell (67420 cells)}

\item{lpjml}{Defines LPJmL version for crop/grass and natveg specific inputs.
Only relevant, if refData = "lpj".}

\item{climatetype}{Switch between different GCM climate scenarios.
Only relevant, if refData = "lpj".}
}
\value{
magpie object in cellular resolution
}
\description{
Calculates the area than can be potentially covered by forests,
             based on environmental conditions.
}
\examples{
\dontrun{
calcOutput("PotentialForestArea", aggregate = FALSE)
}

}
\author{
Patrick v. Jeetze
}
