# Visualizing latent states in a fitted kDGLM model

Visualizing latent states in a fitted kDGLM model

## Usage

``` r
# S3 method for class 'dlm_coef'
plot(
  x,
  var = rownames(x$theta.mean)[x$dynamic],
  cutoff = floor(t/10),
  pred.cred = 0.95,
  plot.pkg = "auto",
  ...
)
```

## Arguments

- x:

  dlm_coef object: The coefficients of a fitted DGLM model.

- var:

  character: The name of the variables to plot (same value passed while
  creating the structure). Any variable whose name partially match this
  variable will be plotted.

- cutoff:

  integer: The number of initial steps that should be skipped in the
  plot. Usually, the model is still learning in the initial steps, so
  the estimated values are not reliable.

- pred.cred:

  numeric: The credibility value for the credibility interval.

- plot.pkg:

  character: A flag indicating if a plot should be produced. Should be
  one of 'auto', 'base', 'ggplot2' or 'plotly'.

- ...:

  Extra arguments passed to the plot method.

## Value

ggplot or plotly object: A plot showing the predictive mean and
credibility interval with the observed data.

## See also

[`fit_model`](fit_model.md),[`coef`](https://rdrr.io/r/stats/coef.html)

Other auxiliary visualization functions for the fitted_dlm class:
[`plot.fitted_dlm()`](plot.fitted_dlm.md),
[`summary.fitted_dlm()`](summary.fitted_dlm.md),
[`summary.searched_dlm()`](summary.searched_dlm.md)

## Examples

``` r
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)

model.coef <- coef(fitted.data)

plot(model.coef)$plot
#> NULL
```
