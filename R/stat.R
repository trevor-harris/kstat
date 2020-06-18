#' Tests if the functions in \code{fmat} and \code{gmat} are equal in distribution
#' @param fmat Matrix of functions. Each column is a function.
#' @param gmat Matrix of functions. Each column is a function. Need to be same length as fmat.
#'
#' @return Value of the K statistic and the associated p value under the null, as a vector.
#'
#' @export
kstat = function(fmat, gmat) {
  ff.xd = int_depth(fmat, fmat)
  fg.xd = int_depth(fmat, gmat)

  gg.xd = int_depth(gmat, gmat)
  gf.xd = int_depth(gmat, fmat)

  ff.cdf = sapply(ff.xd, function(y) mean(ff.xd <= y))
  gf.cdf = sapply(ff.xd, function(y) mean(gf.xd <= y))
  fg.cdf = sapply(gg.xd, function(y) mean(fg.xd <= y))
  gg.cdf = sapply(gg.xd, function(y) mean(gg.xd <= y))

  rate = sqrt((ncol(gmat)*ncol(fmat)) / (ncol(gmat) + ncol(fmat)))

  ksf = max(abs(ff.cdf - gf.cdf))
  ksg = max(abs(gg.cdf - fg.cdf))

  ks = max(ksf, ksg)

  c(ks, 1-ks_cdf(rate*ks))
}


#' Calculates the probability of the value \code{x} under the Kolmogorov Distribution.
#'
#' @param x A positive number.
#' @param n A positive integer. Number of terms to include in the Kolmogorov distribution. Defaults to 20.
#'
#' @return Probability of x.
#'
#' @export
ks_cdf = function(x, n = 20) {
  if(x < 0.05) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}



#' Tests if the functions in \code{fmat} and \code{gmat} are equal in distribution
#' @param fmat Matrix of functions. Each column is a function.
#' @param gmat Matrix of functions. Each column is a function. Need to be same length as fmat.
#' @param perms Positive integer. Number of permutations to construct the approximate permutation distribution.
#'
#'
#' @return Value of the K statistic and the associated p value under the null, as a vector, using a permutation distribution
#'
#' @export
kstat_perm = function(fmat, gmat, perms = 500) {

  # compute KSD statistic
  ks = kstat(fmat, gmat)[1]

  # construct permutation distribution
  hmat = cbind(fmat, gmat)
  hn = ncol(hmat)
  fn = ncol(fmat)

  ksd.dist = rep(0, perms)

  ksd.dist = sapply(1:perms, function(y) {
    hstar = hmat[,sample(1:hn, hn, replace = F)]
    kstat(hstar[,1:fn], hstar[,-(1:fn)])[1]
  })

  # return KS and permutation p value
  c(ks, mean(ksd.dist > ks))
}


