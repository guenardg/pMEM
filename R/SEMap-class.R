## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Spatial eigenvector map class generator and methods**
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
#' SEMap: Spatial EigenMap Class for Predictive Moran's Eigenvector Maps
#' 
#' \code{SEMap} objects store spatial eigenfunctions derived from distance-based
#' weighting matrices, enabling spatially-explicit prediction at sampled and
#' unsampled locations.
#' 
#' @docType class
#' 
#' @name SEMap-class
#' 
#' @aliases SEMap genSEF print.SEMap as.data.frame.SEMap as.matrix.SEMap
#' predict.SEMap
#' 
#' @param x For \code{genSEF}: a numeric matrix or vector of coordinates
#'   (\code{n × d}, where \code{n} is number of sites and \code{d} is dimensions).
#'   For methods: an \code{SEMap}-class object.
#' @param m A distance metric function (e.g., from \code{\link{genDistMetric}}).
#'   Must accept two arguments (\code{x}, \code{y}) and return a numeric or
#'   complex matrix of pairwise distances.
#' @param f A distance weighting function (e.g., from \code{\link{genDWF}}).
#'   Must accept one argument (distances) and return weights of the same type.
#' @param tol Tolerance threshold (positive numeric). Eigenvectors with absolute
#'   eigenvalues below this value are discarded. Default is
#'   \code{.Machine$double.eps^0.5} (~1e-8).
#' @param object An \code{SEMap}-class object (for methods).
#' @param newdata A numeric matrix or vector of coordinates for prediction
#'   (\code{n_new × d}).
#' @param ... Additional arguments passed to methods (e.g., \code{wh} for
#'   selecting specific eigenvectors).
#' @param row.names,optional Passed to \code{as.data.frame}; see
#'   \code{\link[base]{as.data.frame}} for details.
#' 
#' @return 
#' \describe{
#'   \item{\code{genSEF}}{ An \code{SEMap}-class object containing
#'         eigenfunctions, eigenvalues, and prediction methods. }
#'   \item{\code{print.SEMap}}{ \code{NULL} (invisibly); prints summary to
#'         console. }
#'   \item{\code{as.data.frame.SEMap}}{ A \code{data.frame} with eigenvectors as
#'         columns. }
#'   \item{\code{as.matrix.SEMap}}{ A numeric/complex matrix of eigenvectors. }
#'   \item{\code{predict.SEMap}}{ A matrix of eigenfunction scores at
#'         \code{newdata} locations (\code{n_new × k}, where \code{k} is the
#'         number of eigenvectors). }
#' }
#' 
#' @details
#' 
#' Predictive Moran's Eigenvector Maps (pMEM) allows one to model the spatial
#' variability of an environmental variable and use the resulting model for
#' making prediction at any location on and around the sampling points.
#' 
#' \subsection{Algorithm}{
#'   pMEM eigenfunctions are computed as follows:
#'   \enumerate{
#'     \item Compute pairwise distances from coordinates using \code{m(x, x)}.
#'     \item Transform distances to weights using \code{f(d)}.
#'     \item Double-center the weight matrix (row and column means subtracted).
#'     \item Perform eigenvalue decomposition on the centered matrix.
#'     \item Retain eigenvectors with absolute eigenvalues above \code{tol}.
#'   }
#' }
#' 
#' \subsection{Prediction at New Locations}{
#'   The \code{predict} method calculates eigenfunction scores at unsampled
#'   locations by:
#'   \enumerate{
#'     \item Computing distances from training sites to new locations.
#'     \item Transforming to weights and re-centering using training site
#'           centers.
#'     \item Projecting onto the eigenvector basis via matrix multiplication.
#'   }
#' }
#' 
#' \subsection{Complex-Valued Eigenfunctions}{
#'   When using asymmetric distance metrics (via \code{genDistMetric(delta)}),
#'   eigenfunctions are complex-valued. The real and imaginary parts represent
#'   directional spatial patterns (e.g., upstream vs. downstream in rivers).
#' }
#' 
#' \subsection{Modeling Workflow}{
#'   Standard workflow:
#'   \enumerate{
#'     \item Generate \code{SEMap} object with \code{genSEF()}.
#'     \item Extract eigenvectors with \code{as.matrix()} or
#'           \code{as.data.frame()}.
#'     \item Select optimal subset using \code{\link{getMinMSE}} or other
#'           criteria.
#'     \item Fit model (e.g., \code{lm()}, \code{glm()}) with selected
#'           eigenvectors.
#'     \item Predict at new locations using \code{predict(SEMap, newdata)}.
#'   }
#' }
#' 
#' @format An \code{SEMap}-class object is a list with class attribute
#'   \code{"SEMap"} containing the following components:
#' \describe{
#'   \item{\code{show}}{ Function to print object summary. }
#'   \item{\code{getIMoran}}{ Function returning Moran's I coefficients for each
#'         eigenfunction. }
#'   \item{\code{getSEF}}{ Function returning eigenvectors; accepts \code{wh}
#'         argument to select specific columns. }
#'   \item{\code{getLambda}}{ Function returning eigenvalues. }
#'   \item{\code{getPredictor}}{ Function computing eigenfunction scores at new
#'         locations; accepts \code{xx} (coordinates) and \code{wh} (selection). }
#' }
#' 
#' @author \packageAuthor{pMEM}
#' 
#' @importFrom utils head tail
#' 
#' @examples
#' 
#' ## Case 1: one-dimensional symmetrical
#' 
#' n <- 11
#' x <- (n - 1)*seq(0, 1, length.out=n)
#' 
#' ## Store graphical parameters:
#' tmp <- par(no.readonly = TRUE)
#' par(las=1)
#' 
#' sef <- genSEF(x, genDistMetric(), genDWF("Gaussian",3))
#' sef
#' #> A SEMap-class object
#' #> --------------------
#' #> Number of sites: 11
#' #> Directional: no
#' #> Number of components: 10
#' #> Eigenvalues: 3.28700,1.98782,0.79880,...,0.00005,0.00000
#' #> --------------------
#' 
#' ## Extract eigenvectors:
#' dim(as.matrix(sef))  # 11 × 10 matrix
#' #> [1] 11 10
#' 
#' ## Predict at new locations:
#' xx <- (n - 1)*seq(0, 1, 0.01)
#' pred <- predict(sef, xx, wh=1:3)
#' dim(pred)  # 101 × 3 matrix
#' #> [1] 101   3
#' 
#' ## Quick plot of the first eigenfunction:
#' if(interactive()) {
#'   plot(xx, pred[, 1], type="l", xlab="Position", ylab="pMEM_1")
#'   points(x, as.matrix(sef)[, 1], pch=21, bg="black")
#' }
#' 
#' \dontrun{
#' 
#' ## The Second eigenfunction:
#' plot(y = predict(sef, xx, wh=2), x=xx, type="l", ylab="PMEM_2", xlab="x")
#' points(y=as.matrix(sef, wh=2), x=x)
#' 
#' plot(y=predict(sef, xx, wh=5), x=xx, type="l", ylab="PMEM_5", xlab="x")
#' points(y=as.matrix(sef, wh=5), x=x)
#' 
#' ## Case 2: one-dimensional asymmetrical (each has a real and imaginary parts)
#' 
#' sef <- genSEF(x, genDistMetric(delta=pi/8), genDWF("Gaussian",3))
#' 
#' ## First asymmetric eigenfunction:
#' plot(y = Re(predict(sef, xx, wh=1)), x=xx, type="l", ylab="PMEM_1", xlab="x",
#'      ylim=c(-0.35,0.35))
#' lines(y = Im(predict(sef, xx, wh=1)), x=xx, col="red")
#' points(y=Re(as.matrix(sef, wh=1)), x=x)
#' points(y=Im(as.matrix(sef, wh=1)), x=x, col="red")
#' 
#' ## Second asymmetric eigenfunction:
#' plot(y=Re(predict(sef, xx, wh=2)), x=xx, type="l", ylab="PMEM_2", xlab="x",
#'      ylim=c(-0.45,0.35))
#' lines(y = Im(predict(sef, xx, wh=2)), x=xx, col="red")
#' points(y = Re(as.matrix(sef, wh=2)), x=x)
#' points(y = Im(as.matrix(sef, wh=2)), x=x, col="red")
#' 
#' ## Fifth asymmetric eigenfunction:
#' plot(y = Re(predict(sef, xx, wh=5)), x=xx, type="l", ylab="PMEM_5", xlab="x",
#'      ylim=c(-0.45,0.35))
#' lines(y = Im(predict(sef, xx, wh=5)), x=xx, col="red")
#' points(y = Re(as.matrix(sef, wh=5)), x=x)
#' points(y = Im(as.matrix(sef, wh=5)), x=x, col="red")
#' 
#' ## A function to display combinations of the real and imaginary parts:
#' plotAsy <- function(object, xx, wh, a, ylim) {
#'   pp <- predict(object, xx, wh=wh)
#'   plot(y = cos(a)*Re(pp) + sin(a)*Im(pp), x = xx, type = "l",
#'        ylab = "PMEM_5", xlab = "x", ylim=ylim, col="green")
#'   invisible(NULL)
#' }
#' 
#' ## Display combinations at an angle of 45° (pMEM_5):
#' plotAsy(sef, xx, 5, pi/4, ylim=c(-0.45,0.45))
#' 
#' ## Display combinations for other angles:
#' for(i in 0:15)
#'   plotAsy(sef, xx, 5, i*pi/8, ylim=c(-0.45,0.45))
#' 
#' ## Case 3: two-dimensional symmetrical
#' 
#' cbind(
#'   x = c(-0.5,0.5,-1,0,1,-0.5,0.5),
#'   y = c(rep(sqrt(3)/2,2L),rep(0,3L),rep(-sqrt(3)/2,2L))
#' ) -> x2
#' 
#' seq(min(x2[,1L]) - 0.3, max(x2[,1L]) + 0.3, 0.05) -> xx
#' seq(min(x2[,2L]) - 0.3, max(x2[,2L]) + 0.3, 0.05) -> yy
#' 
#' list(
#'   x = xx,
#'   y = yy,
#'   coords = cbind(
#'     x = rep(xx, length(yy)),
#'     y = rep(yy, each = length(xx))
#'   )
#' ) -> ss
#' 
#' cc <- seq(0,1,0.01)
#' cc <- c(rgb(cc,cc,1),rgb(1,1-cc,1-cc))
#' 
#' sef <- genSEF(x2, genDistMetric(), genDWF("Gaussian",3))
#' 
#' scr <- predict(sef, ss$coords)
#' 
#' par(mfrow = c(2,3), mar=0.5*c(1,1,1,1))
#' 
#' for(i in 1L:6) {
#'   image(z=matrix(scr[,i],length(ss$x),length(ss$y)), x=ss$x, y=ss$y, asp=1,
#'         zlim=max(abs(scr[,i]))*c(-1,1), col=cc, axes=FALSE)
#'   points(x = x2[,1L], y = x2[,2L])
#' }
#' 
#' ## Case 4: two-dimensional asymmetrical
#' 
#' sef <- genSEF(x2, genDistMetric(delta=pi/8), genDWF("Gaussian",1))
#' ## Note: default influence angle is 0 (with respect to the abscissa)
#' 
#' ## A function to display combinations of the real and imaginary parts (2D):
#' plotAsy2 <- function(object, ss, a) {
#'   pp <- predict(object, ss$coords)
#'   for(i in 1:6) {
#'     z <- cos(a)*Re(pp[,i]) + sin(a)*Im(pp[,i])
#'     image(z=matrix(z,length(ss$x),length(ss$y)), x=ss$x, y=ss$y, asp=1,
#'           zlim=max(abs(z))*c(-1,1), col=cc, axes=FALSE)
#'   }
#'   invisible(NULL)
#' }
#' 
#' ## Display combinations at an angle of 22°:
#' plotAsy2(sef, ss, pi/8)
#' 
#' ## Combinations at other angles:
#' for(i in 0:23)
#'   plotAsy2(sef, ss, i*pi/12)
#' 
#' ## With an influence of +45° (with respect to the abscissa)
#' sef <- genSEF(x2, genDistMetric(delta=pi/8, theta=pi/4),
#'               genDWF("Gaussian",1))
#' 
#' ## Combinations at other angles:
#' for(i in 0:23)
#'   plotAsy2(sef, ss, i*pi/12)
#' 
#' ## With an influence of +90° (with respect to the abscissa)
#' sef <- genSEF(x2, genDistMetric(delta=pi/8, theta=pi/2),
#'               genDWF("Gaussian",1))
#' 
#' ## Combinations at other angles:
#' for(i in 0:23)
#'   plotAsy2(sef, ss, i*pi/12)
#' 
#' ## With an influence of -45° (with respect to the abscissa)
#' sef <- genSEF(x2, genDistMetric(delta=pi/8, theta=-pi/4),
#'               genDWF("Gaussian",1))
#' 
#' ## Combinations at other angles:
#' for(i in 0:23)
#'   plotAsy2(sef, ss, i*pi/12)
#' 
#' }
#' 
#' ## Reverting to initial graphical parameters:
#' par(tmp)
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib pMEM, .registration = TRUE
#' 
NULL
#' 
#' @describeIn SEMap-class
#' Predictive Moran's Eigenvector Map (pMEM) Generation
#' 
#' Generates a predictive spatial eigenvector map (a SEMap-class object).
#' 
#' @export
genSEF <- function(x, m, f, tol = .Machine$double.eps^0.5) {
  n <- NROW(x)
  d <- m(x,x)
  complex <- is.complex(d)
  w <- f(d)
  if(complex) {
    g <- .Call("pMEM_centerCplx", PACKAGE="pMEM", w, TRUE)
  } else {
    g <- .Call("pMEM_centerReal", PACKAGE="pMEM", w, TRUE)
  }
  eigen(g$centered, symmetric = TRUE) -> eig
  eig$vectors <- eig$vectors[,abs(eig$values) > tol, drop=FALSE]
  eig$values <- eig$values[abs(eig$values) > tol]
  names(eig$values) <- paste("pMEM",1L:length(eig$values),sep="_")
  list(
    if(is.matrix(x)) rownames(x) else names(x),
    names(eig$values)
  ) -> dimnames(eig$vectors)
  p <- eig$vectors %*% diag(eig$values^(-1))
  colnames(p) <- colnames(eig$vectors)
  structure(
    list(
      show = function() {
        cat("A SEMap-class object\n--------------------\n")
        cat(sprintf("Number of sites: %d\n",n))
        cat(sprintf("Directional: %s\n",if(complex) "yes" else "no"))
        cat(sprintf("Number of components: %d\n",length(eig$values)))
        ev <- sprintf("%.5f", eig$values)
        if(length(eig$values) > 6L)
          ev <- c(head(ev,3L),"...",tail(ev,2))
        cat(sprintf("Eigenvalues: %s\n", paste(ev, collapse=",")))
        cat("--------------------\n")
        invisible(NULL)
      },
      getIMoran = function() {
        I <- (eig$values - 1)*n/(sum(w) - n)
        if(complex) Re(I) else I
      },
      getSEF = function(wh)
        if(missing(wh)) eig$vectors else eig$vectors[,wh,drop=FALSE],
      getLambda = function() eig$values,
      getPredictor = function(xx, wh) {
        dd <- m(x,xx)
        ww <- f(dd)
        if(complex) {
          gg <- .Call("pMEM_recenterCplx", PACKAGE="pMEM", ww, g$centers, TRUE)
        } else {
          gg <- .Call("pMEM_recenterReal", PACKAGE="pMEM", ww, g$centers, TRUE)
        }
        if(missing(wh)) gg %*% p else gg %*% p[,wh,drop=FALSE]
      }
    ),
    class = "SEMap"
  )
}
#'
#'@describeIn SEMap-class
#' Print SEMap-class
#' 
#' A print method for \code{SEMap-class} objects.
#' 
#' @method print SEMap
#' @export
print.SEMap <- function(x, ...)
  x$show()
#' 
#' @describeIn SEMap-class
#' An \code{as.data.frame} Method for \code{SEMap-class} Objects
#' 
#' A method to extract the spatial eigenvectors from an \code{SEMap-class}
#' object as a data frame.
#' 
#' @method as.data.frame SEMap
#' @export
as.data.frame.SEMap <- function(x, row.names = NULL, optional = FALSE, ...)
  as.data.frame(x$getSEF(...), row.names = row.names, optional = optional)
#' 
#' @describeIn SEMap-class
#' An \code{as.matrix} Method for \code{SEMap-class} Objects
#' 
#' A method to extract the spatial eigenvectors from an \code{SEMap-class}
#' object as a matrix.
#' 
#' @method as.matrix SEMap
#' @export
as.matrix.SEMap <- function(x, ...)
  x$getSEF(...)
#' 
#' @describeIn SEMap-class
#' A \code{predict} Method for \code{SEMap-class} Objects
#' 
#' A method to obtain predictions from an \code{SEMap-class} object.
#' 
#' @method predict SEMap
#' @export
predict.SEMap <- function(object, newdata, ...)
  object$getPredictor(newdata, ...)
#' 
