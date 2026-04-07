# coef.fitted_dlm

Evaluates the predictive values for the observed values used to fit the
model and its latent states. Predictions can be made with smoothed
values, with filtered values or h-steps ahead.

## Usage

``` r
# S3 method for class 'fitted_dlm'
coef(
  object,
  t.eval = seq_len(object$t),
  lag = -1,
  pred.cred = 0.95,
  eval.pred = FALSE,
  eval.metric = FALSE,
  ...
)
```

## Arguments

- object:

  fitted_dlm: The fitted model to be use for evaluation.

- t.eval:

  numeric: A vector of positive integers indicating the time index from
  which to extract predictions. The default is to extract to evaluate
  the model at all observed times.

- lag:

  integer: The relative offset for forecast. Values for time t will be
  calculated based on the filtered values of time t-h. If lag is
  negative, then the smoothed distribution for the latent states will be
  used.

- pred.cred:

  numeric: The credibility level for the C.I..

- eval.pred:

  boolean: A flag indicating if the predictions should be calculated.

- eval.metric:

  boolean: A flag indicating if the model density (f(M\|y)) should be
  calculated. Only used when lag\<0.

- ...:

  Extra arguments passed to the coef method.

## Value

A list containing:

- data data.frame: A table with the model evaluated at each observed
  time.

- theta.mean matrix: The mean of the latent states at each time.
  Dimensions are n x t, where t is the size of t.eval and n is the
  number of latent states.

- theta.cov array: A 3D-array containing the covariance matrix of the
  latent states at each time. Dimensions are n x n x t, where t is the
  size of t.eval and n is the number of latent states.

- lambda.mean matrix: The mean of the linear predictor at each time.
  Dimensions are k x t, where t is the size of t.eval and k is the
  number of linear predictors.

- lambda.cov array: A 3D-array containing the covariance matrix for the
  linear predictor at each time. Dimensions are k x k x t, where t is
  the size of t.eval and k is the number of linear predictors.

- log.like, mae, mase, rae, mse, interval.score: The metric value at
  each time.

- conj.param list: A list containing, for each outcome, a data.frame
  with the parameter of the conjugated distribution at each time.

## See also

Other auxiliary functions for fitted_dlm objects:
[`eval_dlm_norm_const()`](https://silvaneojunior.github.io/kDGLM/reference/eval_dlm_norm_const.md),
[`fit_model()`](https://silvaneojunior.github.io/kDGLM/reference/fit_model.md),
[`forecast.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/forecast.fitted_dlm.md),
[`kdglm()`](https://silvaneojunior.github.io/kDGLM/reference/kdglm.md),
[`simulate.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/simulate.fitted_dlm.md),
[`smoothing()`](https://silvaneojunior.github.io/kDGLM/reference/smoothing.md),
[`update.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/update.fitted_dlm.md)

## Examples

``` r
# Poisson case
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data = data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)

var.vals <- coef(fitted.data)
```
