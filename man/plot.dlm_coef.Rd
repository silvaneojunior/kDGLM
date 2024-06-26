% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_helper.R
\name{plot.dlm_coef}
\alias{plot.dlm_coef}
\title{Visualizing latent states in a fitted kDGLM model}
\usage{
\method{plot}{dlm_coef}(
  x,
  var = rownames(x$theta.mean)[x$dynamic],
  cutoff = floor(t/10),
  pred.cred = 0.95,
  plot.pkg = "auto",
  ...
)
}
\arguments{
\item{x}{dlm_coef object: The coefficients of a fitted DGLM model.}

\item{var}{character: The name of the variables to plot (same value passed while creating the structure). Any variable whose name partially match this variable will be plotted.}

\item{cutoff}{integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.}

\item{pred.cred}{numeric: The credibility value for the credibility interval.}

\item{plot.pkg}{character: A flag indicating if a plot should be produced. Should be one of 'auto', 'base', 'ggplot2' or 'plotly'.}

\item{...}{Extra arguments passed to the plot method.}
}
\value{
ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
}
\description{
Visualizing latent states in a fitted kDGLM model
}
\examples{

data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)

model.coef <- coef(fitted.data)

plot(model.coef)$plot

}
\seealso{
\code{\link{fit_model}},\code{\link{coef}}

Other auxiliary visualization functions for the fitted_dlm class: 
\code{\link{plot.fitted_dlm}()},
\code{\link{summary.fitted_dlm}()},
\code{\link{summary.searched_dlm}()}
}
\concept{auxiliary visualization functions for the fitted_dlm class}
