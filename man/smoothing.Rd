% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{smoothing}
\alias{smoothing}
\title{Auxiliary function for model smoothing}
\usage{
smoothing(model)
}
\arguments{
\item{model}{A fitted_dlm object.}
}
\value{
A fitted_dlm object with smoothed means (mts) and covariance matrix (Cts) for each observation.
}
\description{
Auxiliary function for model smoothing
}
\seealso{
Other auxiliary functions for fitted_dlm objects: 
\code{\link{coef.fitted_dlm}()},
\code{\link{eval_dlm_norm_const}()},
\code{\link{fit_model}()},
\code{\link{forecast.fitted_dlm}()},
\code{\link{kdglm}()},
\code{\link{simulate.fitted_dlm}()},
\code{\link{update.fitted_dlm}()}
}
\concept{auxiliary functions for fitted_dlm objects}
