% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_gamma.R
\name{convert_Normal_Gamma}
\alias{convert_Normal_Gamma}
\title{convert_Normal_Gamma}
\usage{
convert_Normal_Gamma(conj.param, parms)
}
\arguments{
\item{conj.param}{list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).}

\item{parms}{list: A list of extra known parameters of the distribution. Not used in this function.}
}
\value{
The parameters of the Normal distribution of the linear predictor.
}
\description{
Calculates the parameters of the log-Normal that best approximates the given Inverse-Gamma distribution.
The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
}
\seealso{
Other auxiliary functions for a Gamma outcome with known shape: 
\code{\link{convert_Gamma_Normal}()},
\code{\link{gamma_pred}()},
\code{\link{update_Gamma}()}
}
\concept{auxiliary functions for a Gamma outcome with known shape}
\keyword{internal}
