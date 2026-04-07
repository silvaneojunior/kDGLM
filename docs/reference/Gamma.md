# Gamma outcome for kDGLM models

Creates an outcome with gamma distribution with the chosen parameters
(can only specify 2).

## Usage

``` r
Gamma(
  phi = NA,
  mu = NA,
  alpha = NA,
  beta = NA,
  sigma = NA,
  data,
  offset = as.matrix(data)^0
)
```

## Arguments

- phi:

  character or numeric: Name of the linear predictor associated with the
  shape parameter of the gamma distribution. If numeric, this parameter
  is treated as known and equal to the value passed. If a character, the
  parameter is treated as unknown and equal to the exponential of the
  associated linear predictor. It cannot be specified with alpha.

- mu:

  character: Name of the linear predictor associated with the mean
  parameter of the gamma distribution. The parameter is treated as
  unknown and equal to the exponential of the associated linear
  predictor.

- alpha:

  character: Name of the linear predictor associated with the shape
  parameter of the gamma distribution. The parameter is treated as
  unknown and equal to the exponential of the associated linear
  predictor. It cannot be specified with phi.

- beta:

  character: Name of the linear predictor associated with the rate
  parameter of the gamma distribution. The parameter is treated as
  unknown and equal to the exponential of the associated linear
  predictor. It cannot be specified with sigma.

- sigma:

  character: Name of the linear predictor associated with the scale
  parameter of the gamma distribution. The parameter is treated as
  unknown and equal to the exponential of the associated linear
  predictor. It cannot be specified with beta.

- data:

  numeric: Values of the observed data.

- offset:

  numeric: The offset at each observation. Must have the same shape as
  data.

## Value

An object of the class dlm_distr

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
[`Multinom()`](https://silvaneojunior.github.io/kDGLM/reference/Multinom.md),
[`Normal()`](https://silvaneojunior.github.io/kDGLM/reference/Normal.md),
[`Poisson()`](https://silvaneojunior.github.io/kDGLM/reference/Poisson.md),
[`summary.dlm_distr()`](https://silvaneojunior.github.io/kDGLM/reference/summary.dlm_distr.md)

## Examples

``` r
structure <- polynomial_block(mu = 1, D = 0.95)

Y <- (cornWheat$corn.log.return[1:500] - mean(cornWheat$corn.log.return[1:500]))**2
outcome <- Gamma(phi = 0.5, mu = "mu", data = Y)
fitted.data <- fit_model(structure, corn = outcome)
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     corn: Gamma
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : 3508.409
#> Interval Score        : 0.00197
#> Mean Abs. Scaled Error: 0.93721
#> ---
plot(fitted.data, plot.pkg = "base")

```
