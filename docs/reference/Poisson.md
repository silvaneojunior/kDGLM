# Poisson outcome for kDGLM models

Creates an outcome with Poisson distribution with the chosen parameter.

## Usage

``` r
Poisson(lambda, data, offset = as.matrix(data)^0)
```

## Arguments

- lambda:

  character: The name of the linear predictor associated with the rate
  (mean) parameter of the Poisson distribution. The parameter is treated
  as unknown and equal to the exponential of the associated linear
  predictor.

- data:

  numeric: The values of the observed data.

- offset:

  numeric: The offset at each observation. Must have the same shape as
  data.

## Value

A object of the class dlm_distr

## Details

For evaluating the posterior parameters, we use the method proposed in
Alves et al. (2024) .

For the details about the implementation see dos Santos et al. (2024) .

## References

Mariane Branco Alves, Helio S. Migon, Raíra Marotta, Junior, Silvaneo
Vieira dos Santos (2024). “k-parametric Dynamic Generalized Linear
Models: a sequential approach via Information Geometry.” 2201.05387.  
  
Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”

## See also

[`fit_model`](https://silvaneojunior.github.io/kDGLM/reference/fit_model.md)

Other auxiliary functions for a creating outcomes:
[`Gamma()`](https://silvaneojunior.github.io/kDGLM/reference/Gamma.md),
[`Multinom()`](https://silvaneojunior.github.io/kDGLM/reference/Multinom.md),
[`Normal()`](https://silvaneojunior.github.io/kDGLM/reference/Normal.md),
[`summary.dlm_distr()`](https://silvaneojunior.github.io/kDGLM/reference/summary.dlm_distr.md)

## Examples

``` r
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, D = 0.95, order = 2)
season <- harmonic_block(rate = 1, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data = data)

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
#> Log-likelihood        : -625.2975
#> Interval Score        : 138.36923
#> Mean Abs. Scaled Error:   0.68377
#> ---

plot(fitted.data, plot.pkg = "base")

```
