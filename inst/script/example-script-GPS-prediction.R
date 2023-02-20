## Example script for GPS prediction: 1 parameter -> planes in space (G_{2,3}).
library(gpsr)
library(Matrix)

## Problem setup ------------------------------------------------------------
## 2-dim subspace in 3-dim Euclidean space.
n <- 3
k <- 2

## Initial subspace (at t = 0)
normal0 <- c(1, 1, 1) #Initial normal vector
Q <- svd(matrix(normal0, ncol = 1), nu = n)$u
X0 <- Q[,-1] #Initial orthonormal basis

#' Rotation maxtrix of a given angle around z-axis.
zRotate <- function(theta) {
    c <- cos(theta)
    s <- sin(theta)
    matrix(c(c, s, 0, -s, c, 0, 0, 0, 1), nrow = n)
}

#' True subspace-valued function.
#' @param t time, the input parameter.
#' @param ometa Rotation speed in rad / time, defaults to pi/2.
#' @return Orthonormal basis of the rotated subspace.
subspaceFn <- function(t, omega = pi/2) {
    zRotate(omega * t) %*% X0
}

## Training data
l <- 4L #Number of training points
thetaTrain <- round(seq(0, 100, length.out = l)) / 100
listXtrain <- lapply(thetaTrain, subspaceFn)

## Test data
nTest <- 120
thetaTest <- seq(0, 1, length.out = nTest)
listXtest <- lapply(thetaTest, subspaceFn)

## GPS prediction ------------------------------------------------------------
## (0) Prepare X
AX <- gpsr:::listMatrix2Array(listXtrain)
str(AX)
X <- Matrix::Matrix(AX, n, k * l) #Concatenated bases.

## (1) Pre-processing
retGPSprep <- gpsr::GPSubspacePrepSVD(X)
list2env(retGPSprep, .GlobalEnv) #Put variables in global environment.

## (2) Hyperparameter tuning
## Require (thetaTrain, XtX, VbtX; Jk)
Jk <- gpsr:::J(k)
## Option 1: Default length-scale, no tuning.
(len <- gpsr::defaultLength(d = 1, l))
## Option 2: Use LOOCV error.
lenUpper <- len * 2
lenLower <- len * 0.5
system.time(ret <- optimize(gpsr::hSSDist, lower = lenLower, upper = lenUpper, tol = 1e-2))
## Option 3: Use LOOCV error and gradient.
ToleranceLevel <- function(x) list(factr = 10^(15-x))
tol4 <- ToleranceLevel(4)
system.time(ret2 <- optim(len, gpsr::hSSDist, gpsr::gSSDist, method = "L-BFGS-B",
                          lower = lenLower, upper = lenUpper, control = tol4))
## Not optimal implementation: `optim()` needs separate arguments for objective and gradient,
## while `gSSDist()` can compute both.

## (3) Prediction
len <- ret$minimum
## Correlation matrix
K <- gpsr::kerSEmat(thetaTrain, len = len)
system.time(ret <- purrr::map(thetaTest, ~gpsr::GPSubspacePredEVD(., thetaTrain, len, K)))
## (Optional) explicit form of reduced bases.
listVpred <- purrr::map(ret, ~svdX$u %*% .$Vcirc)

## Error analysis --------------------------------------------------

## Error measures: subspace distance.
## Distance between true function and GPS predictions.
distGPS <- purrr::map2_dbl(listVpred, listXtest, gpsr:::RiemannianDist)
plot(thetaTest, distGPS, type = 'o', main = "GPS prediction error in (normalized) subspace distance.")
abline(h = 0, lty = 2)
points(thetaTrain, rep(0, l), pch = 19, col = "blue") #Training points.

## Representative standardized normal vectors.
#' Standardized normal vector of a plane in space.
getNormal <- function(X) svd(X, nu = n)$u[,n]
#' Sign of inner product
signInnerProd <- function(x, y) sign(sum(x * y))
#' Representative normal vectors.
#' @param Mn Matrix of normal vectors.
getStandardNormals <- function(Mn) {
    posSign <- purrr::map_dbl(seq(ncol(Mn) - 1), ~signInnerProd(Mn[,.], Mn[,.+1]))
    posSign0 <- signInnerProd(normal0, Mn[,1])
    flipSign <- cumprod(c(posSign0, posSign))
    return(Mn %*% diag(flipSign))
}
#' Representative standardized normal vectors from a list of basis.
getStdNormals <- function(lsX){
    Mn <- vapply(lsX, getNormal, double(n))
    getStandardNormals(Mn)
}

stdNormalTrain <- getStdNormals(listXtrain)
stdNormalTest <- getStdNormals(listXtest)
stdNormalPred <- getStdNormals(listVpred)

dev.new()
matplot(thetaTest, t(stdNormalTest), type = 'l', lty = 1)
matplot(thetaTest, t(stdNormalPred), type = 'l', lty = 2, add = TRUE)
matpoints(thetaTrain, t(stdNormalTrain), pch = 1)
title("Coordinates of normal vectors")
txt <- c(paste(c('x', 'y', 'z'), "true (training points)"), paste(c('x', 'y', 'z'), "predict"))
legend("bottomleft", legend = txt, lty = rep(c(1, 2), each = 3),
       col = rep(seq(3), 2), pch = rep(c(1, NA), each = 3))
