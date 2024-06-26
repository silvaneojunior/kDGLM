% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_helper.R
\name{if.na}
\alias{if.na}
\title{if.na}
\usage{
if.na(vec, val)
}
\arguments{
\item{vec}{A vector or matrix.}

\item{val}{The value to replace NA with.}
}
\value{
A vector or matrix with the same dimensions as the input, where any NA values have been replaced by the specified val argument.
}
\description{
This function is wrapper for ifelse(is.na(vec),vec,val)
}
\keyword{internal}
