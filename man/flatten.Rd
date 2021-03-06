% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{flatten}
\alias{flatten}
\title{Flatten a 3D array into a matrix by converting each matrix in the first two indicies into a vector.
Mostly used to feed spatial data into the kstat function.}
\usage{
flatten(arr)
}
\arguments{
\item{arr}{Three dimensional array.}
}
\value{
Two dimensional representation of \code{arr}
}
\description{
Flatten a 3D array into a matrix by converting each matrix in the first two indicies into a vector.
Mostly used to feed spatial data into the kstat function.
}
