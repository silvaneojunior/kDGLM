# coefficients.fitted_dlm

This method is wrapper for the coef method.

## Usage

``` r
# S3 method for class 'fitted_dlm'
coefficients(object, ...)
```

## Arguments

- object:

  A fitted_dlm object.

- ...:

  Arguments passed to coef.

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

[`coef.fitted_dlm`](https://silvaneojunior.github.io/kDGLM/reference/coef.fitted_dlm.md)
