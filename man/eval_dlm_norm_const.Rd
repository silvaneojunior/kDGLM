% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{eval_dlm_norm_const}
\alias{eval_dlm_norm_const}
\title{Auxiliary function for evaluating normalizing constant for the posterior of a fitted DLM.}
\usage{
eval_dlm_norm_const(model, lin.pred = FALSE, safe.mode)
}
\arguments{
\item{model}{fitted_dlm: A fitted_dlm object.}

\item{lin.pred}{boolean: A flag indicating if the normalizing constant should be calculated using the linear predictors.}
}
\value{
A scalar representing the normalizing constant for the posterior of a fitted DLM.
}
\description{
Evaluates the normalizing constant for the posterior of a fitted DLM.
}
\examples{

data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data = data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)
eval_dlm_norm_const(fitted.data)

}
\seealso{
Other auxiliary functions for fitted_dlm objects: 
\code{\link{coef.fitted_dlm}()},
\code{\link{fit_model}()},
\code{\link{forecast.fitted_dlm}()},
\code{\link{kdglm}()},
\code{\link{simulate.fitted_dlm}()},
\code{\link{smoothing}()},
\code{\link{update.fitted_dlm}()}
}
\concept{auxiliary functions for fitted_dlm objects}
\keyword{internal}
