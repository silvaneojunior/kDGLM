# Zero sum prior

Defines the prior of a structural block to be such that the latent
states sum zero with probability one.

## Usage

``` r
zero_sum_prior(
  block,
  var.index = 1:block$n,
  weights = rep(1, length(var.index))
)
```

## Arguments

- block:

  dlm_block object: The structural block.

- var.index:

  integer: The index of the variables from which to set the prior.

- weights:

  numeric: A vector indicating which linear transformation of the data
  is 0 with probability 1. Default is equivalent to a zero-sum
  restriction.

## Value

A dlm_block object with the desired prior.

## Details

The covariance matrix of the evolution and the drift parameter are also
altered to guarantee that the zero sum condition will always hold. The
discount factor must be the same for all variables whose prior is being
modified. For the details about the implementation see dos Santos et al.
(2024) .

## References

Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”

## See also

Other auxiliary functions for defining priors.:
[`CAR_prior()`](CAR_prior.md), [`joint_prior()`](joint_prior.md)

## Examples

``` r
polynomial_block(mu = 1, D = 0.95) |>
  block_mult(5) |>
  zero_sum_prior()
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
