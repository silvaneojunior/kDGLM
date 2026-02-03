# Summary for a fitted kDGLM model

Prints a report for a fitted_dlm object.

## Usage

``` r
# S3 method for class 'fitted_dlm'
summary(
  object,
  t = object$t,
  lag = -1,
  metric.lag = 1,
  metric.cutoff = floor(object$t/10),
  pred.cred = 0.95,
  ...
)
```

## Arguments

- object:

  A fitted_dlm object.

- t:

  Integer: The time index for the latent states.

- lag:

  Integer: The number of steps ahead used for the evaluating the latent
  states. Use lag\<0 for the smoothed distribution, If lag==0 for the
  filtered distribution and lag=h for the h-step-ahead prediction.

- metric.lag:

  Integer: The number of steps ahead used for the evaluating the
  predictions used when calculating metrics. Use metric.lag\<0 for the
  smoothed distribution, If metric.lag==0 for the filtered distribution
  and metric.lag=h for the h-step-ahead prediction.

- metric.cutoff:

  Integer: The cutoff time index for the metric calculation. Values
  before that time will be ignored.

- pred.cred:

  numeric: The credibility interval to be used for the interval score.

- ...:

  Extra arguments passed to the coef method.#'

## Value

No return value, called to print a summary of the fitted kDGLM model.

## See also

Other auxiliary visualization functions for the fitted_dlm class:
[`plot.dlm_coef()`](plot.dlm_coef.md),
[`plot.fitted_dlm()`](plot.fitted_dlm.md),
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
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     AirPassengers: Poisson
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : -580.2514
#> Interval Score        : 129.48462
#> Mean Abs. Scaled Error:   0.47013
#> ---
```
