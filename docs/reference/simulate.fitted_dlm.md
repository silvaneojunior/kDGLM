# Draw samples from the distribution of the latent states

This is function draws samples from the latent states using the backward
sampling algorithm. See West and Harrison (1997) , chapter 15, for
details.

## Usage

``` r
# S3 method for class 'fitted_dlm'
simulate(object, nsim, seed = NULL, lag = -1, ...)
```

## Arguments

- object:

  fitted_dlm: A fitted model from which to sample.

- nsim:

  integer: The number of samples to draw.

- seed:

  integer: An object specifying if and how the random number generator
  should be initialized.

- lag:

  integer: The relative offset for forecast. Values for time t will be
  calculated based on the filtered values of time t-h. If lag is
  negative, then the smoothed distribution for the latent states will be
  used.

- ...:

  Extra arguments passed to the plot method.

## Value

A list containing the following values:

- theta array: An array containing a sample of the latent states.
  Dimensions are n x t x nsim, where n is the number of latent states in
  the model and t is the number of observed values.

- lambda array: An array containing a sample of the linear predictors.
  Dimensions are k x t x nsim, where k is the number of linear
  predictors in the model and t is the number of observed values.

- param list: A named list containing, for each model outcome, an array
  with the samples of the parameters of the observational model. Each
  array will have dimensions l x t x nsim, where l is the number of
  parameters in the observational model and t is the number of observed
  values.

## See also

Other auxiliary functions for fitted_dlm objects:
[`coef.fitted_dlm()`](coef.fitted_dlm.md),
[`eval_dlm_norm_const()`](eval_dlm_norm_const.md),
[`fit_model()`](fit_model.md),
[`forecast.fitted_dlm()`](forecast.fitted_dlm.md),
[`kdglm()`](kdglm.md), [`smoothing()`](smoothing.md),
[`update.fitted_dlm()`](update.fitted_dlm.md)

## Examples

``` r
structure <- polynomial_block(mu = 1, D = 0.95) +
  polynomial_block(V = 1, D = 0.95)

outcome <- Normal(mu = "mu", V = "V", data = cornWheat$corn.log.return[1:500])
fitted.data <- fit_model(structure, corn = outcome)

sample <- simulate(fitted.data, 5000)
```
