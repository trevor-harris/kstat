#' Tests if the functions in \code{gmat} have the same distribution as the functions in \code{fmat}
#'
#' @param gmat Matrix of functions. Each column is a function.
#' @param fmat Matrix of functions. Each column is a function. Need to be same length as gmat.
#'
#' @return Value of the K statistic and the associated p value under the null, as a vector.
#'
#' @export
kstat = function(fmat, gmat) {
  ff.xd = xdepth(f, f)
  fg.xd = xdepth(f, g)

  gg.xd = xdepth(g, g)
  gf.xd = xdepth(g, f)

  ff.cdf = sapply(ff.xd, function(y) mean(ff.xd <= y))
  gf.cdf = sapply(ff.xd, function(y) mean(gf.xd <= y))
  fg.cdf = sapply(gg.xd, function(y) mean(fg.xd <= y))
  gg.cdf = sapply(gg.xd, function(y) mean(gg.xd <= y))

  rate = sqrt((ncol(g)*ncol(f)) / (ncol(g) + ncol(f)))

  ksf = max(abs(ff.cdf - gf.cdf))
  ksg = max(abs(gg.cdf - fg.cdf))

  ks = max(ksf, ksg)

  c(ks, 1-ks_cdf(rate*ks))
}


#' Calculates the probability of the value \code{x} under the Kolmogorov Distribution.
#'
#' @param x A positive number.
#' @param n A positive number. Number of terms to include in the Kolmogorov distribution. Defaults to 20.
#'
#' @return Probability of x.
#'
#' @export
ks_cdf = function(x, n = 20) {
  if(x < 0.05) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}
