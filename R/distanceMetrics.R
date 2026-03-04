## **************************************************************************
##
##    (c) 2023-2026 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Distance Metric Generator Function**
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
##    R source code file
##
## **************************************************************************
##
#' Generate a Distance Metric Function
#'
#' \code{genDistMetric} returns a function that calculates pairwise distances
#' between rows of coordinate matrices, with optional asymmetry parameters.
#' 
#' @param delta Optional asymmetry parameter. When provided, the returned
#'   distance metric becomes complex-valued, with modulus equal to the Euclidean
#'   distance and argument determined by \code{delta}. For 1D data, the argument
#'   is ±\code{delta}; for 2D data, it is the cosine of the angular difference
#'   relative to \code{theta}.
#' @param theta Influence angle (in radians) for 2D asymmetric metrics. Used
#'   only when \code{delta} is provided. Default is \code{0}.
#' 
#' @returns A function with signature \code{function(x, y)} that returns a
#' numeric matrix of pairwise distances. When \code{y} is omitted, the function
#' returns a square symmetric (or Hermitian) matrix of distances among rows of
#' \code{x}. For asymmetric metrics (\code{delta} provided), the matrix is
#' Hermitian with strictly real eigenvalues.
#' 
#' @details
#' \subsection{Symmetric vs. Asymmetric Metrics}{
#'   When \code{delta} is missing (default), the returned function computes the
#'   standard Euclidean distance. When \code{delta} is provided, the metric
#'   becomes complex-valued: its modulus equals the Euclidean distance, and its
#'   argument is determined by \code{delta} as follows:
#'   \itemize{
#'     \item \strong{1D data (transects):} The argument is \code{+delta} for
#'           pairs where the second point is downstream (higher coordinate) and
#'           \code{-delta} otherwise.
#'     \item \strong{2D data:} The argument is the cosine of the angular
#'           difference between the bearing from point A to B and the influence
#'           angle \code{theta}.
#'   }
#'   In both cases, the distance from A → B has the opposite argument sign as
#'   B → A, ensuring the pairwise distance matrix is Hermitian with strictly
#'   real eigenvalues.
#' }
#' 
#' \subsection{Function Generator Pattern}{
#'   \code{genDistMetric} uses a function-generator pattern: it returns a
#'   closure that embeds the specified \code{delta} and \code{theta} values
#'   in its environment. To change these parameters, call
#'   \code{genDistMetric} again with new arguments.
#' }
#' 
#' @author \packageAuthor{pMEM}
#' 
#' @examples
#' ## A five-point equidistant transect:
#' n <- 5
#' x <- (n - 1)*seq(0, 1, length.out=n)
#' 
#' ## Symmetric (Euclidean) metric function:
#' mSym <- genDistMetric()
#' 
#' ## The pairwise symmetric metric between the rows of x:
#' mSym(x)
#' #>      [,1] [,2] [,3] [,4] [,5]
#' #> [1,]    0    1    2    3    4
#' #> [2,]    1    0    1    2    3
#' #> [3,]    2    1    0    1    2
#' #> [4,]    3    2    1    0    1
#' #> [5,]    4    3    2    1    0
#' 
#' ## Distances between x and a finer grid xx:
#' xx <- (n - 1)*seq(0, 1, 0.05)
#' mSym(x, xx)  # Rectangular matrix: 5 rows × 21 columns
#' 
#' ## Asymmetric metric with delta = 0.2:
#' mAsy <- genDistMetric(delta = 0.2)
#' mAsy(x)  # Hermitian matrix (complex values)
#' Mod(mAsy(x))  # Modulus equals Euclidean distances
#' 
#' ## Asymmetric distances between x and xx:
#' mAsy(x, xx)
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib pMEM, .registration = TRUE
#' 
#' @export
genDistMetric <- function(delta, theta = 0)
  if(missing(delta)) {
    function(x, y) {
      if(!is.matrix(x)) x <- as.matrix(x)
      if(missing(y)) y <- x
      if(!is.matrix(y)) y <- as.matrix(y)
      .Call("pMEM_EuclidReal", PACKAGE="pMEM", x, y)
    }
  } else {
    function(x, y) {
      if(!is.matrix(x)) x <- as.matrix(x)
      if(missing(y)) y <- x
      if(!is.matrix(y)) y <- as.matrix(y)
      .Call("pMEM_EuclidCplx2D", PACKAGE="pMEM", x, y, delta, theta)
    }
  }
#'
