# This file is part of MSTWeatherGen package.
#
# MSTWeatherGen is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# MSTWeatherGen is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with MSTWeatherGen.
# If not, see <https://www.gnu.org/licenses/>.


#' @encoding UTF-8
#' @title Multivariate Space-Time Weather Generator
#' @description Simulation model of multivariate spatio-temporal meteorological variables.
#' It provides tools for estimating the model parameters and generating
#` synthetic weather data that can be used for a variety of applications,
#` including climate research, agricultural or hydrological modeling.
#' @aliases MSTWeatherGen-package MSTWeatherGen
#'
#' @author Said Obakrim \email{saidobak@@gmail.com}
#' @author Jean-François Rey \email{jean-francois@@inrae.fr}
#'
#' Maintainer: Said Obakrim \email{saidobak@@gmail.com}
#' @docType package
#' @name MSTWeatherGen-package
#' @details \tabular{ll}{
#'          Package: \tab MSTWeatherGen\cr
#'          Type: \tab Package\cr
#'          Version: \tab 0.1.0\cr
#'          Date: \tab 2024-03-30\cr
#'          License: \tab GPL (>=3)\cr
#'          }
#'
#' @keywords internal
#' @references
#' ## When referencing the simulation model, please cite the following article:
#'
#'
#' ## When referencing the R package, please cite the following package:
#'
#' url: https://cran.r-project.org/package=MSTWeatherGen.
#' @examples
#' \dontrun{
#' library("MSTWeather")
#'
#' }

#' @import stats
#' @importFrom Matrix nearPD
#' @import lubridate
#' @import ggplot2
#' @importFrom dplyr group_by tally mutate case_when summarise
#' @importFrom VGAM pbinorm
#' @import crch
#' @import abind
#' @import viridis
#' @import stringr
#' @import PTAk
#' @import mclust
#' @importFrom utils packageVersion
#'
"_PACKAGE"

# @title Print package information
# @name getInfo
# @description Displays some information about the package
# @importFrom utils packageVersion
getInfo <- function() {
    packageStartupMessage("Package: MSTWeatherGen | MST Weather Generator")
    packageStartupMessage("Version: ", appendLF = FALSE)
    packageStartupMessage(utils::packageVersion("MSTWeatherGen"))
    packageStartupMessage("License: GPL (>= 3)")
}

# @title Things to do at package attach
# @name .onAttach
# @param libname a character string giving the library directory where
#  the package defining the namespace was found.
# @param pkgname a character string giving the name of the package.
# @description Print package information and check dependencies
.onAttach <- function(libname, pkgname) {
    getInfo()
}

