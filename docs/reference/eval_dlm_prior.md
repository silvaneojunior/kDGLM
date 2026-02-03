# Auxiliary function for evaluating the prior density of a DLM

Evaluates the prior density for a set of parameters theta in a DLM. The
structure of the DLM is taken to be that of the fitted_dlm object passed
as input.

## Usage

``` r
eval_dlm_prior(theta, model, lin.pred = FALSE)
```

## Arguments

- theta:

  matrix: A matrix representing the set of parameter for which to
  evaluate the density. Its size should be n x t, where n is the number
  of latent states and t is the length of the time series;

- model:

  fitted_dlm object: A fitted_dlm object.

- lin.pred:

  boolean: A flag indicating if theta represents the linear predictors.

## Value

A scalar representing the log density evaluated at theta.

## Examples

``` r
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data = data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)
eval_dlm_prior(fitted.data$mts, fitted.data)
#> [1] 3628.202
```
