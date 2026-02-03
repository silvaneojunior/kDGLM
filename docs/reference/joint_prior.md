# Joint prior

Defines the joint prior of a structural block.

## Usage

``` r
joint_prior(
  block,
  var.index = 1:block$n,
  a1 = block$a1[var.index],
  R1 = block$R1[var.index, var.index]
)
```

## Arguments

- block:

  dlm_block object: The structural block.

- var.index:

  Integer: The index of the variables from which to set the prior.

- a1:

  Numeric: The prior mean.

- R1:

  Matrix: The prior covariance matrix.

## Value

A dlm_block object with the desired prior.

## Details

The discount factor must be the same for all variables whose prior is
being modified. For the details about the implementation see dos Santos
et al. (2024) .

## References

Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”

## See also

Other auxiliary functions for defining priors.:
[`CAR_prior()`](CAR_prior.md), [`zero_sum_prior()`](zero_sum_prior.md)

## Examples

``` r
polynomial_block(mu = 1, D = 0.95) |>
  block_mult(5) |>
  joint_prior(var.index = 1:2, R1 = matrix(c(1, 0.5, 0.5, 1), 2, 2))
#> Mixed DLM block.
#> latent states: 
#>     Var.Poly.1: Level (1 variable(s))
#>     Var.Poly.2: Level (1 variable(s))
#>     Var.Poly.3: Level (1 variable(s))
#>     Var.Poly.4: Level (1 variable(s))
#>     Var.Poly.5: Level (1 variable(s))
#> 
#> Linear predictors: 
#>     mu.1
#>     mu.2
#>     mu.3
#>     mu.4
#>     mu.5
#> 
#> Status: defined
#> Serie length: 1
#> Interventions at: 
#> Number of latent states: 5
#> Number of linear predictors: 5
```
