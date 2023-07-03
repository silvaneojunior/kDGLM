% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_normal.R
\name{Normal}
\alias{Normal}
\title{Normal outcome for kDGLM models}
\usage{
Normal(
  mu,
  Tau = NA,
  Sigma = NA,
  Sd = NA,
  outcome,
  offset = outcome^0,
  alt_method = FALSE
)
}
\arguments{
\item{mu}{character: Name of the linear predictor associated with the mean parameter of the Normal distribution. The parameter is treated as unknowed and equal to the associated linear predictor.}

\item{Tau}{character or numeric: If Tau is a character, it is interpreted as the names of the linear preditors associated with the precisions parameter of the Normal distribution. If Tau is numeric, the precision is considered known and equal to the value of Tau, otherwise, the precison is considere unknowned and equal to the exponential of the linear predictor informed in Tau. If the outcome is a Multivariate Normal, then Tau must be a matrix and, if the precision is unknowned, the elements outside it's main diagonal are treated as the linear predictor associated with the correlation between the each coordinate of the outcome, otherwise Tau is treated as the precision matrix. The user cannot be specify Tau with Sigma or Sd.}

\item{Sigma}{character or numeric: If Sigma is a character, it is interpreted as the names of the linear preditors associated with the variance parameter of the Normal distribution. If Sigma is numeric, the variance is considered known and equal to the value of Sigma, otherwise, the variance is considere unknowned and equal to the exponential of the linear predictor informed in Sigma. If the outcome is a Multivariate Normal, then Sigma must be a matrix and, if the variance is unknowned, the elements outside it's main diagonal are treated as the linear predictor associated with the correlation between the each coordinate of the outcome, otherwise Sigma is treated as the covariance matrix. The user cannot be specify Sigma with Tau or Sd.}

\item{Sd}{character or numeric: If Sd is a character, it is interpreted as the names of the linear preditors associated with the standard deviation parameter of the Normal distribution. If Sd is numeric, the standard deviation is considered known and equal to the value of Sd, otherwise, the precison is considere unknowned and equal to the exponential of the linear predictor informed by in Sd. If the outcome is a Multivariate Normal, then Tau must be a matrix and the elements outside it's main diagonal are treated as the correlation (or the name of the linear predictor associated) between the each coordinate of the outcome. The user cannot be specify Sd with Tau or Sigma.}

\item{outcome}{vector: Values of the observed data.}

\item{offset}{vector: The offset at each observation. Must have the same shape as outcome.}
}
\value{
A object of the class dlm_distr
}
\description{
Creates an outcome with Normal distribution with the chosen parameters (can only specify 2,).
}
\details{
For evaluating the posterior parameters, we use the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.

For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
}
\examples{

# Normal case
T <- 200
mu <- rnorm(T, 0, 0.1)
data <- rnorm(T, cumsum(mu))

level <- polynomial_block(
  mu = 1,
  D = 0.95
)
variance <- polynomial_block(
  sigma2 = 1
)

# Known variance
outcome <- Normal(mu = "mu", Sigma = 1, outcome = data)

fitted_data <- fit_model(level, outcomes = outcome)
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot

# Unknown variance
outcome <- Normal(mu = "mu", Sigma = "sigma2", outcome = data)

fitted_data <- fit_model(level, variance, outcomes = outcome)
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{fit_model}}

Other {auxiliary functions for a creating outcomes}: 
\code{\link{Gamma}()},
\code{\link{Multinom}()},
\code{\link{Poisson}()}
}
\concept{{auxiliary functions for a creating outcomes}}