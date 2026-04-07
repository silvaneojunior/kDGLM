# update.fitted_dlm

update.fitted_dlm

## Usage

``` r
# S3 method for class 'fitted_dlm'
update(object, ...)
```

## Arguments

- object:

  fitted_dlm: The fitted model to be updated.

- ...:

  Extra variables necessary for updating (covariates, observed values,
  etc.).

## Value

A fitted_dlm object.

## Details

If an a covariate is necessary for updating, it should be passed as a
named argument. Its name must follow this structure: \<block
name\>.Covariate\<.index\>. If there is only one pulse in the associated
block the index is omitted. If an a pulse is necessary for updating, it
should be passed as a named argument. Its name must follow this
structure: \<block name\>.Pulse\<.index\>. If there is only one pulse in
the associated block the index is omitted. If an offset is necessary for
updating, it should be passed along with the observed data. See example.

## See also

Other auxiliary functions for fitted_dlm objects:
[`coef.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/coef.fitted_dlm.md),
[`eval_dlm_norm_const()`](https://silvaneojunior.github.io/kDGLM/reference/eval_dlm_norm_const.md),
[`fit_model()`](https://silvaneojunior.github.io/kDGLM/reference/fit_model.md),
[`forecast.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/forecast.fitted_dlm.md),
[`kdglm()`](https://silvaneojunior.github.io/kDGLM/reference/kdglm.md),
[`simulate.fitted_dlm()`](https://silvaneojunior.github.io/kDGLM/reference/simulate.fitted_dlm.md),
[`smoothing()`](https://silvaneojunior.github.io/kDGLM/reference/smoothing.md)

## Examples

``` r
level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

# Only first 100 observations (for the sake of the example)
outcome <- Poisson(lambda = "rate", data = c(AirPassengers)[1:100])

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)

updated.fit <- update(fitted.data, AirPassengers = list(data = c(AirPassengers)[101:144]))
# If a offset was present, the user should pass its value when updating
# updated.fit=update(fitted.data,
#                     AirPassengers=list(
#                      data=c(AirPassengers)[101:144],
#                      offset= ... ))
```
