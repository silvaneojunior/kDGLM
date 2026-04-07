# Auxiliary function for evaluating normalizing constant for the posterior of a fitted DLM.

Evaluates the normalizing constant for the posterior of a fitted DLM.

## Usage

``` r
eval_dlm_norm_const(model, lin.pred = FALSE)
```

## Arguments

- model:

  fitted_dlm: A fitted_dlm object.

- lin.pred:

  boolean: A flag indicating if the normalizing constant should be
  calculated using the linear predictors.

## Value

A scalar representing the normalizing constant for the posterior of a
fitted DLM.

## See also

Other auxiliary functions for fitted_dlm objects:
[`coef.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/coef.fitted_dlm.md),
[`fit_model()`](https://silvaneojunior.github.io/kDGLM/reference/fit_model.md),
[`forecast.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/forecast.fitted_dlm.md),
[`kdglm()`](https://silvaneojunior.github.io/kDGLM/reference/kdglm.md),
[`simulate.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/simulate.fitted_dlm.md),
[`smoothing()`](https://silvaneojunior.github.io/kDGLM/reference/smoothing.md),
[`update.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/update.fitted_dlm.md)

## Examples

``` r
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data = data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)
eval_dlm_norm_const(fitted.data)
#> [1] -598.9172
```
