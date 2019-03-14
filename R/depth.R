#' Compute the univariate tukey depth of a function \code{g} with respect to a set of functions \code{fmat}
#'
#' @param g Vector. Function that you want to find the univariate tukey depths of.
#' @param fmat Matrix of functions. Each column is a function. Need to be same length as g.
#'
#' @return Vector of pointwise tukey depths
#'
#' @export
depth = function(g, fmat) {

  fn = ncol(fmat)
  depth = rep(0, length(g))

  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }

  return(depth)
}

#' Compute the integrated tukey depth of each function in  \code{gmat} with respect to a set of functions \code{fmat}
#'
#' @param gmat Matrix of functions. Each column is a function.
#' @param fmat Matrix of functions. Each column is a function. Need to be same length as gmat.
#'
#' @return Vector of integrated tukey depths
#'
int_depth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}
