## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **St. Marguerite river salmon data set**
##
##    This file is part of pMEM
##
##    pMEM is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    pMEM is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with pMEM. If not, see <https://www.gnu.org/licenses/>.
##
##    Data set documentation generation file
##
## **************************************************************************
##
#' Sainte-Marguerite River Atlantic Salmon Parr Data
#' 
#' Spatially-referenced observations of juvenile Atlantic salmon (parr) density
#' and habitat variables along a 1520 m transect in the Sainte-Marguerite River,
#' Québec, Canada.
#' 
#' @docType data
#' 
#' @name salmon
#' 
#' @format A \code{data.frame} with 76 rows (sampling sites) and 5 columns:
#' \describe{
#'   \item{Position}{Numeric: distance along transect from \code{Bardsville}
#'     (m).}
#'   \item{Abundance}{Integer: count of Atlantic salmon parr observed per 20 m
#'     site.}
#'   \item{Depth}{Numeric: mean water depth (m) at thalweg.}
#'   \item{Velocity}{Numeric: mean current velocity (m/s) at thalweg.}
#'   \item{Substrate}{Numeric: mean substrate grain size (mm, D50) at thalweg.}
#' }
#' 
#' @details
#' \subsection{Sampling Design}{
#'   Data were collected on July 7, 2002, along a 1520 m river segment starting
#'   at \code{Bardsville} (48°23'01.59'' N, 70°12'10.05'' W). The transect was
#'   divided into 76 contiguous 20 m sites. At each site, snorkelers recorded
#'   parr abundance while moving upstream in a zigzag pattern to minimize
#'   disturbance.
#' }
#' 
#' \subsection{Habitat Measurements}{
#'   Environmental variables were measured at the thalweg (channel centerline)
#'   at the midpoint of each 20 m site (i.e., 10 m from each boundary):
#'   \itemize{
#'     \item \strong{Depth}: Mean water depth (m)
#'     \item \strong{Velocity}: Mean current velocity (m/s)
#'     \item \strong{Substrate}: Mean grain size (mm, D50 metric)
#'   }
#' }
#' 
#' \subsection{Coordinate System}{
#'   \code{Position} is a local Cartesian coordinate (meters) increasing
#'   upstream from the \code{Bardsville} reference point. No projected CRS is
#'   assigned; treat as a 1D transect for spatial modeling.
#' }
#' 
#' @source Daniel Boisclair, Département de sciences biologiques, Université de
#' Montréal, Montréal, Québec, Canada.
#' 
#' @references
#' Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010. Multiscale
#' codependence analysis: an integrated approach to analyse relationships across
#' scales. Ecology 91: 2952-2964 <doi:10.1890/09-0460.1>
#' 
#' @seealso
#' Bouchard, J. and Boisclair, D. 2008. The relative importance of local,
#' lateral, and longitudinal variables on the development of habitat quality
#' models for a river. Can. J. Fish. Aquat. Sci. 65: 61-73 <doi:10.1139/f07-140>
#' 
#' @examples
#' data(salmon)
#' 
#' ## Quick summary:
#' summary(salmon)
#' #>     Position        Abundance            Depth         Velocity       Substrate
#' #> Min.   :1280   Min.   : 0.000   Min.   :0.3600   Min.   :0.1900   Min.   : 51.0
#' #> 1st Qu.:1655   1st Qu.: 0.000   1st Qu.:0.5075   1st Qu.:0.4700   1st Qu.:159.5
#' #> Median :2030   Median : 2.000   Median :0.6050   Median :0.5750   Median :220.0
#' #> Mean   :2030   Mean   : 2.092   Mean   :0.6184   Mean   :0.5868   Mean   :213.9
#' #> 3rd Qu.:2405   3rd Qu.: 3.000   3rd Qu.:0.6900   3rd Qu.:0.6725   3rd Qu.:265.2
#' #> Max.   :2780   Max.   :14.000   Max.   :1.0600   Max.   :1.0700   Max.   :381.0
#' 
#' ## Simple plot of parr abundance along the transect:
#' if(interactive()) {
#'   plot(salmon$Position, salmon$Abundance, type = "b", pch = 21, bg = "black",
#'        xlab = "Position along transect (m)", ylab = "Parr abundance")
#' }
#' 
#' @keywords salmon spatial river transect habitat parr
#' 
NULL
