# CAR prior

Defines the prior of a structural block as a Conditional Autoregressive
(CAR) prior.

## Usage

``` r
CAR_prior(
  block,
  adj.matrix,
  scale,
  rho,
  sum.zero = FALSE,
  var.index = 1:block$n
)
```

## Arguments

- block:

  dlm_block object: The structural block.

- adj.matrix:

  matrix: The adjacency matrix.

- scale:

  numeric: The tau parameter for the CAR model (see references).

- rho:

  numeric: The rho parameter for the CAR model (see references).

- sum.zero:

  Bool: If true, all latent states will add to 0.

- var.index:

  integer: The index of the variables from which to set the prior.

## Value

A dlm_block object with the desired prior.

## Details

The filtering algorithm used in this package requires a proper prior for
the latent space. As such, this implementation of the CAR prior imposes
a zero-sum constraint in the regional effects. The discount factor must
be the same for all variables whose prior is being modified.

For a revision of the CAR prior, see Schmidt and Nobre (2018) .

For the details about the implementation see dos Santos et al. (2024) .

## References

Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”  
  
Alexandra M. Schmidt, Widemberg S. Nobre (2018). “Conditional
Autoregressive (CAR) Model.” In *Wiley StatsRef: Statistics Reference
Online*, chapter Conditional Autoregressive (CAR) Model, 1-11. John
Wiley & Sons, Ltd. ISBN 9781118445112,
[doi:10.1002/9781118445112.stat08048](https://doi.org/10.1002/9781118445112.stat08048)
,
https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118445112.stat08048,
<https://onlinelibrary.wiley.com/doi/abs/10.1002/9781118445112.stat08048>.

## See also

Auxiliary functions for creating structural blocks
[`polynomial_block`](https://silvaneojunior.github.io/kDGLM/reference/polynomial_block.md),
[`regression_block`](https://silvaneojunior.github.io/kDGLM/reference/regression_block.md),
[`harmonic_block`](https://silvaneojunior.github.io/kDGLM/reference/harmonic_block.md),
[`TF_block`](https://silvaneojunior.github.io/kDGLM/reference/tf_block.md).

Other auxiliary functions for defining priors.:
[`joint_prior()`](https://silvaneojunior.github.io/kDGLM/reference/joint_prior.md),
[`zero_sum_prior()`](https://silvaneojunior.github.io/kDGLM/reference/zero_sum_prior.md)

## Examples

``` r
# Creating an arbitrary adjacency matrix
adj.matrix <- matrix(
  c(
    0, 1, 1, 0, 0,
    1, 0, 1, 0, 0,
    1, 1, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 1, 0
  ),
  5, 5,
  byrow = TRUE
)

polynomial_block(mu = 1, D = 0.95) |>
  block_mult(5) |>
  CAR_prior(scale = 9, rho = 1, adj.matrix = adj.matrix)
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
