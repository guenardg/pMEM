## **************************************************************************
##
##    (c) 2023-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Simple orthogonal term selection regression**
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
#' Orthogonal Term Selection via Cross-Validated MSE
#' 
#' \code{getMinMSE} performs forward selection of orthogonal predictors by
#' minimizing out-of-sample mean squared error (MSE) on a testing dataset.
#' 
#' @param U Training matrix of orthogonal predictors (\code{n_train × p}),
#'   where \code{n_train} is the number of training observations and \code{p}
#'   is the number of predictors (e.g., spatial eigenvectors). Must be
#'   orthonormal (columns have unit norm and are mutually orthogonal).
#' @param y Training response vector (length \code{n_train}).
#' @param Up Testing matrix of orthogonal predictors (\code{n_test × p}),
#'   evaluated at testing locations.
#' @param yy Testing response vector (length \code{n_test}).
#' @param complete If \code{TRUE} (default), return full selection results
#'   (all MSE values, coefficient order, etc.). If \code{FALSE}, return only
#'   the minimum MSE and associated coefficient threshold.
#' 
#' @returns 
#' If \code{complete = TRUE}, a list with the following elements:
#' \describe{
#'   \item{betasq}{ Squared standardized regression coefficients (length
#'   \code{p}). }
#'   \item{nullmse}{ Null model MSE: variance of testing responses around their
#'   mean. }
#'   \item{mse}{ MSE values for incremental models (length \code{p + 1}); 
#'   \code{mse[1]} is the null MSE. }
#'   \item{ord}{ Indices of predictors sorted by decreasing absolute coefficient
#'   value. }
#'   \item{wh}{ Index of the best model (1-based). Value \code{1} means the null 
#'   model is best; value \code{k + 1} means the first \code{k} predictors from
#'   \code{ord} are selected. }
#' }
#' 
#' If \code{complete = FALSE}, a list with two elements:
#' \describe{
#'   \item{betasq}{ Squared standardized coefficient at the optimal model. }
#'   \item{mse}{ Minimum MSE value achieved. }
#' }
#' 
#' @details
#' This function allows one to calculate a simple model, involving only
#' the spatial eigenvectors and a single response variable. The coefficients
#' are estimated on a training data set; the ones that are retained are chosen
#' on the basis of minimizing the mean squared error on the testing data set.
#' 
#' \subsection{Algorithm}{
#'   The selection procedure proceeds as follows:
#'   \enumerate{
#'     \item Compute regression coefficients: \code{b = t(U) \%*\% y} (valid
#'           because \code{U} is orthonormal).
#'     \item Sort coefficients by decreasing absolute value; store order in
#'           \code{ord}.
#'     \item Compute null MSE: variance of \code{yy} around \code{mean(y)}.
#'     \item For each predictor (in sorted order), add its partial prediction to
#'           the model and compute out-of-sample MSE on \code{yy}.
#'     \item Identify the model with minimum MSE; return corresponding
#'           coefficient threshold and MSE.
#'   }
#' }
#' 
#' \subsection{Orthonormality Requirement}{
#'   Predictors in \code{U} must be orthonormal (columns have unit norm and are 
#'   mutually orthogonal). This condition is met by design for spatial
#'   eigenvectors produced by \code{genSEF()}, but users must ensure it for
#'   custom predictors.
#' }
#' 
#' \subsection{Performance}{
#'   This function is implemented in C++ via Rcpp for efficiency. For
#'   \code{p = 50} predictors and \code{n = 100} observations, it is
#'   approximately 50× faster than a pure R implementation.
#' }
#' 
#' @author \packageAuthor{pMEM}
#' 
#' @examples
#' ## Loading the 'salmon' dataset
#' data("salmon")
#' seq(1,nrow(salmon),3) -> test      # Indices of the testing set.
#' (1:nrow(salmon))[-test] -> train   # Indices of the training set.
#' 
#' ## A set of locations located 1 m apart:
#' xx <- seq(min(salmon$Position) - 20, max(salmon$Position) + 20, 1)
#' 
#' ## Lists to contain the results:
#' mseRes <- list()
#' sel <- list()
#' lm <- list()
#' prd <- list()
#' 
#' ## Generate the spatial eigenfunctions:
#' genSEF(
#'   x = salmon$Position[train],
#'   m = genDistMetric(),
#'   f = genDWF("Gaussian",40)
#' ) -> sefTrain
#' 
#' ## Spatially-explicit modelling of the channel depth:
#' 
#' ## Calculate the minimum MSE model:
#' getMinMSE(
#'   U = as.matrix(sefTrain),
#'   y = salmon$Depth[train],
#'   Up = predict(sefTrain, salmon$Position[test]),
#'   yy = salmon$Depth[test]
#' ) -> mseRes[["Depth"]]
#' 
#' ## This is the coefficient of prediction:
#' 1 - mseRes$Depth$mse[mseRes$Depth$wh]/mseRes$Depth$nullmse
#' 
#' ## Storing graphical parameters:
#' tmp <- par(no.readonly = TRUE)
#' 
#' ## Changing the graphical margins:
#' par(mar=c(4,4,2,2))
#' 
#' ## Plot of the MSE values:
#' plot(mseRes$Depth$mse, type="l", ylab="MSE", xlab="order", axes=FALSE,
#'      ylim=c(0.005,0.025))
#' points(x=1:length(mseRes$Depth$mse), y=mseRes$Depth$mse, pch=21, bg="black")
#' axis(1)
#' axis(2, las=1)
#' abline(h=mseRes$Depth$nullmse, lty=3)  # Dotted line: the null MSE
#' 
#' ## A list of the selected spatial eigenfunctions:
#' sel[["Depth"]] <- sort(mseRes$Depth$ord[1:mseRes$Depth$wh])
#' 
#' ## A linear model built using the selected spatial eigenfunctions:
#' lm(
#'   formula = y~.,
#'   data = cbind(
#'     y = salmon$Depth[train],
#'     as.data.frame(sefTrain, wh=sel$Depth)
#'   )
#' ) -> lm[["Depth"]]
#' 
#' ## Calculating predictions of depth at each 1 m intervals:
#' predict(
#'   lm$Depth,
#'   newdata = as.data.frame(
#'     predict(
#'       object = sefTrain,
#'       newdata = xx,
#'       wh = sel$Depth
#'     )
#'   )
#' ) -> prd[["Depth"]]
#' 
#' ## Plot of the predicted depth (solid line), and observed depth for the
#' ## training set (black markers) and testing set (red markers):
#' plot(x=xx, y=prd$Depth, type="l", ylim=range(salmon$Depth, prd$Depth), las=1,
#'      ylab="Depth (m)", xlab="Location along the transect (m)")
#' points(x = salmon$Position[train], y = salmon$Depth[train], pch=21,
#'        bg="black")
#' points(x = salmon$Position[test], y = salmon$Depth[test], pch=21, bg="red")
#' 
#' ## Prediction of the velocity, substrate, and using them to predict the parr
#' ## density.
#' \dontrun{
#' 
#' ## Spatially-explicit modelling of the water velocity:
#' 
#' ## Calculate the minimum MSE model:
#' getMinMSE(
#'   U = as.matrix(sefTrain),
#'   y = salmon$Velocity[train],
#'   Up = predict(sefTrain, salmon$Position[test]),
#'   yy = salmon$Velocity[test]
#' ) -> mseRes[["Velocity"]]
#' 
#' ## This is the coefficient of prediction:
#' 1 - mseRes$Velocity$mse[mseRes$Velocity$wh]/mseRes$Velocity$nullmse
#' 
#' ## Plot of the MSE values:
#' plot(mseRes$Velocity$mse, type="l", ylab="MSE", xlab="order", axes=FALSE,
#'      ylim=c(0.010,0.030))
#' points(x=1:length(mseRes$Velocity$mse), y=mseRes$Velocity$mse, pch=21,
#'        bg="black")
#' axis(1)
#' axis(2, las=1)
#' abline(h=mseRes$Velocity$nullmse, lty=3)
#' 
#' ## A list of the selected spatial eigenfunctions:
#' sel[["Velocity"]] <- sort(mseRes$Velocity$ord[1:mseRes$Velocity$wh])
#' 
#' ## A linear model build using the selected spatial eigenfunctions:
#' lm(
#'   formula = y~.,
#'   data = cbind(
#'     y = salmon$Velocity[train],
#'     as.data.frame(sefTrain, wh=sel$Velocity)
#'   )
#' ) -> lm[["Velocity"]]
#' 
#' ## Calculating predictions of velocity at each 1 m intervals:
#' predict(
#'   lm$Velocity,
#'   newdata = as.data.frame(
#'     predict(
#'       object = sefTrain,
#'       newdata = xx,
#'       wh = sel$Velocity
#'     )
#'   )
#' ) -> prd[["Velocity"]]
#' 
#' ## Plot of the predicted velocity (solid line), and observed velocity for the
#' ## training set (black markers) and testing set (red markers):
#' plot(x=xx, y=prd$Velocity, type="l",
#'      ylim=range(salmon$Velocity, prd$Velocity),
#'      las=1, ylab="Velocity (m/s)", xlab="Location along the transect (m)")
#' points(x = salmon$Position[train], y = salmon$Velocity[train], pch=21,
#'        bg="black")
#' points(x = salmon$Position[test], y = salmon$Velocity[test], pch=21,
#'        bg="red")
#' 
#' ## Spatially-explicit modelling of the mean substrate size (D50):
#' 
#' ## Calculate the minimum MSE model:
#' getMinMSE(
#'   U = as.matrix(sefTrain),
#'   y = salmon$Substrate[train],
#'   Up = predict(sefTrain, salmon$Position[test]),
#'   yy = salmon$Substrate[test]
#' ) -> mseRes[["Substrate"]]
#' 
#' ## This is the coefficient of prediction:
#' 1 - mseRes$Substrate$mse[mseRes$Substrate$wh]/mseRes$Substrate$nullmse
#' 
#' ## Plot of the MSE values:
#' plot(mseRes$Substrate$mse, type="l", ylab="MSE", xlab="order", axes=FALSE,
#'      ylim=c(1000,6000))
#' points(x=1:length(mseRes$Substrate$mse), y=mseRes$Substrate$mse, pch=21,
#'        bg="black")
#' axis(1)
#' axis(2, las=1)
#' abline(h=mseRes$Substrate$nullmse, lty=3)
#' 
#' ## A list of the selected spatial eigenfunctions:
#' sel[["Substrate"]] <- sort(mseRes$Substrate$ord[1:mseRes$Substrate$wh])
#' 
#' ## A linear model build using the selected spatial eigenfunctions:
#' lm(
#'   formula = y~.,
#'   data = cbind(
#'     y = salmon$Substrate[train],
#'     as.data.frame(sefTrain, wh=sel$Substrate)
#'   )
#' ) -> lm[["Substrate"]]
#' 
#' ## Calculating predictions of D50 at each 1 m intervals:
#' predict(
#'   lm$Substrate,
#'   newdata = as.data.frame(
#'     predict(
#'       object = sefTrain,
#'       newdata = xx,
#'       wh = sel$Substrate
#'     )
#'   )
#' ) -> prd[["Substrate"]]
#' 
#' ## Plot of the predicted D50 (solid line), and observed D50 for the training
#' ## set (black markers) and testing set (red markers):
#' plot(x=xx, y=prd$Substrate, type="l",
#'      ylim=range(salmon$Substrate, prd$Substrate), las=1, ylab="D50 (mm)",
#'      xlab="Location along the transect (m)")
#' points(x = salmon$Position[train], y = salmon$Substrate[train], pch=21,
#'        bg="black")
#' points(x = salmon$Position[test], y = salmon$Substrate[test], pch=21,
#'        bg="red")
#' 
#' ## Spatially-explicit modelling of Atlantic salmon parr abundance using
#' ## x=channel depth + water velocity + D50 + pMEM:
#' 
#' ## Requires suggested package glmnet to perform elasticnet regression:
#' 
#' library(glmnet)
#' 
#' ## Calculation of the elastic net model (cross-validated):
#' cv.glmnet(
#'   y = salmon$Abundance[train],
#'   x = cbind(
#'     Depth = salmon$Depth[train],
#'     Velocity = salmon$Velocity[train],
#'     Substrate = salmon$Substrate[train],
#'     as.matrix(sefTrain)
#'   ),
#'   family = "poisson"
#' ) -> cvglm
#' 
#' ## Calculating predictions for the test data:
#' predict(
#'   cvglm,
#'   newx = cbind(
#'     Depth = salmon$Depth[test],
#'     Velocity = salmon$Velocity[test],
#'     Substrate = salmon$Substrate[test],
#'     predict(sefTrain, salmon$Position[test])
#'   ),
#'   s="lambda.min",
#'   type = "response"
#' ) -> yhatTest
#' 
#' ## Calculating predictions for the transect (1 m separated data):
#' predict(
#'   cvglm,
#'   newx = cbind(
#'     Depth = prd$Depth,
#'     Velocity = prd$Velocity,
#'     Substrate = prd$Substrate,
#'     predict(sefTrain, xx)
#'   ),
#'   s = "lambda.min",
#'   type = "response"
#' ) -> yhatTransect
#' 
#' ## Plot of the predicted Atlantic salmon parr abundance (solid line, with the
#' ## depth, velocity, and D50 also predicted using spatially-explicit
#' ## submodels), the observed abundances for the training set (black markers),
#' ## the observed abundances for the testing set (red markers), and the
#' ## predicted abundances for the testing set calculated on the basis of
#' ## observed depth, velocity, and median substrate grain size:
#' plot(x=xx, y=yhatTransect, type="l",
#'      ylim=range(salmon$Abundance,yhatTransect), las=1,
#'      ylab="Abundance (fish)", xlab="Location along the transect (m)")
#' points(x=salmon$Position[train], y=salmon$Abundance[train], pch=21,
#'        bg="black")
#' points(x=salmon$Position[test], y=salmon$Abundance[test], pch=21, bg="red")
#' points(x=salmon$Position[test], y=yhatTest, pch=21, bg="green")
#' 
#' }
#' ## Restoring previous graphical parameters:
#' par(tmp)
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib pMEM, .registration = TRUE
#' 
#' @export
getMinMSE <- function(U, y, Up, yy, complete = TRUE) {
  if(is.complex(U)) {
    stop("Not yet implemented for complex numbers!")
  } else
    .Call("pMEM_getMinMSEReal", PACKAGE="pMEM", U, y, Up, yy, complete)
}
#'
