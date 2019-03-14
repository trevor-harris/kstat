#' Generates realizations of one dimensional Gaussian process with exponential covariance function
#'
#' @param n Positive integer. Number of functions to generate.
#' @param mu numeric value or vector of length \code{pts} Mean value of the Gaussian process. Defaults to 0.
#' @param sd Positive numeric or vector of length \code{pts} Standard deviation of the Gaussian process. Defaults to 1.
#' @param r Positive numeric or vector of length \code{pts} Range of the Gaussian process. Defaults to 10.
#' @param pts Positive integer. Grid size to sample each functions on. Defaults to 50.
#'
#' @note Setting a single number in \code{mu}, \code{sd}, or \code{r} means that the parameter will be constant over all observation
#' points. Setting it to a vector allows you to vary it by observation site. Vector length must be equal to \code{pts}.
#'
#' @return Matrix of functions. Each column is a separate function.
#'
#' @export
gp1d = function(n = 100, mu = 0, sd = 1, r = 10, pts = 50) {
  grid = 1:pts
  distmat = as.matrix(dist(grid))

  # calc sigma with cov kernel
  sigma = exp(-distmat / r)

  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)

  gps = matrix(0, pts, n)
  for(f in 1:n) {
    gps[,f] = (sigma.half %*% rnorm(pts, sd = sd)) + mu
  }
  return(gps)
}

#' Generates realizations of two dimensional Gaussian process with exponential covariance function
#'
#' @param n Positive integer. Number of functions to generate.
#' @param mu numeric value or matrix of size \code{pts^2} Mean value of the Gaussian process. Defaults to 0.
#' @param sd Positive numeric or vector of length \code{pts^2} Standard deviation of the Gaussian process. Defaults to 1.
#' @param r Positive numeric or vector of length \code{pts^2} Range of the Gaussian process. Defaults to 30.
#' @param pts Positive integer. Functions will be sampled on a \code{pts} x \code{pts} grid. Defaults to 30.
#'
#' @note Setting a single number in \code{mu}, \code{sd}, or \code{r} means that the parameter will be constant over all observation
#' points. Setting it to a vector allows you to vary it by observation site. Vector length must be equal to \code{pts}.
#'
#' @return 3D array of functions. First two indicies correspond to values in the 2D GP. Third index is the function number.
#'
#' @export
gp2d = function(n = 100, mu = 0, sd = 1, l = 30, pts = 30) {
  grid = 1:pts
  grid = expand.grid(grid, grid)
  distmat = as.matrix(dist(grid))

  # calc sigma with cov kernel
  sigma = exp(-distmat / l)

  sigma.eig = eigen(sigma)
  sigma.half = sigma.eig$vectors %*% diag(sqrt(sigma.eig$values)) %*% t(sigma.eig$vectors)

  gps = array(0, dim=c(pts, pts, fields))
  for(f in 1:fields) {
    gps[,,f] = (sigma.half %*% rnorm(pts^2, sd = sd)) + mu
  }
  return(gps)
}

#' Flatten a 3D array into a matrix by converting each matrix in the first two indicies into a vector.
#' Mostly used to feed spatial data into the kstat function.
#'
#' @param arr Three dimensional array.
#'
#' @return Two dimensional representation of \code{arr}
#'
#' @export
flatten = function(arr) {
  matrix(arr, prod(dim(arr)[1:2]), dim(arr)[3])
}
