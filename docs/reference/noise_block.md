# noise_block

Creates the structure for a Noise block. This block represents an
independent random noise that should be added to the linear predictor.
The variance of the noise cannot be formally estimated, as such we use a
discount strategy similar to that of West and Harrison (1997) to specify
it.

## Usage

``` r
noise_block(..., name = "Noise", D = 0.99, R1 = 0.1, H = 0)

noise(name = "Noise", D = 0.99, R1 = 0.1, H = 0, X = 1)
```

## Arguments

- ...:

  Named values for the planning matrix.

- name:

  String: An optional argument providing the name for this block. Can be
  useful to identify the models with meaningful labels, also, the name
  used will be used in some auxiliary functions.

- D:

  scalar or vector: A sequence of values specifying the desired discount
  factor for each time. It should have length 1 or t, where t is the
  size of the series. If both D and H are specified, the value of D is
  ignored.

- R1:

  scalar: The prior variance of the noise.

- H:

  scalar: The variance of the noise. If both D and H are specified, the
  value of D is ignored.

- X:

  Vector or scalar: An argument providing the values of the covariate
  X_t.

## Value

A dlm_block object containing the following values:

- FF Array: A 3D-array containing the regression matrix for each time.
  Its dimension should be n x k x t, where n is the number of latent
  states, k is the number of linear predictors in the model and t is the
  time series length.

- FF.labs Matrix: A n x k character matrix describing the type of value
  of each element of FF.

- G Matrix: A 3D-array containing the evolution matrix for each time.
  Its dimension should be n x n x t, where n is the number of latent
  states and t is the time series length.

- G.labs Matrix: A n x n character matrix describing the type of value
  of each element of G.

- D Array: A 3D-array containing the discount factor matrix for each
  time. Its dimension should be n x n x t, where n is the number of
  latent states and t is the time series length.

- H Array: A 3D-array containing the covariance matrix of the noise for
  each time. Its dimension should be the same as D.

- a1 Vector: The prior mean for the latent vector.

- R1 Matrix: The prior covariance matrix for the latent vector.

- var.names list: A list containing the variables indexes by their name.

- order Positive integer: Same as argument.

- n Positive integer: The number of latent states associated with this
  block (2).

- t Positive integer: The number of time steps associated with this
  block. If 1, the block is compatible with blocks of any time length,
  but if t is greater than 1, this block can only be used with blocks of
  the same time length.

- k Positive integer: The number of outcomes associated with this block.
  This block can only be used with blocks with the same outcome length.

- pred.names Vector: The name of the linear predictors associated with
  this block.

- monitoring Vector: The combination of monitoring, monitoring and
  monitoring.pulse.

- type Character: The type of block (Noise).

## Details

For the details about the implementation see dos Santos et al. (2024) .

For the details about dynamic regression models in the context of DLMs,
see West and Harrison (1997) , chapters 6 and 9.

## References

Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”  
  
Mike West, Jeff Harrison (1997). *Bayesian Forecasting and Dynamic
Models (Springer Series in Statistics)*. Springer-Verlag. ISBN
0387947256.

## See also

[`fit_model`](https://silvaneojunior.github.io/kDGLM/reference/fit_model.md)

Other auxiliary functions for structural blocks:
[`TF_block()`](https://silvaneojunior.github.io/kDGLM/reference/tf_block.md),
[`block_mult()`](https://silvaneojunior.github.io/kDGLM/reference/block_mult.md),
[`block_rename()`](https://silvaneojunior.github.io/kDGLM/reference/block_rename.md),
[`block_superpos()`](https://silvaneojunior.github.io/kDGLM/reference/block_superpos.md),
[`ffs_block()`](https://silvaneojunior.github.io/kDGLM/reference/ffs_block.md),
[`harmonic_block()`](https://silvaneojunior.github.io/kDGLM/reference/harmonic_block.md),
[`intervention()`](https://silvaneojunior.github.io/kDGLM/reference/intervention.md),
[`polynomial_block()`](https://silvaneojunior.github.io/kDGLM/reference/polynomial_block.md),
[`regression_block()`](https://silvaneojunior.github.io/kDGLM/reference/regression_block.md),
[`specify.dlm_block()`](https://silvaneojunior.github.io/kDGLM/reference/specify.dlm_block.md),
[`summary.dlm_block()`](https://silvaneojunior.github.io/kDGLM/reference/summary.dlm_block.md)

## Examples

``` r
noise_block(mu = 1, D = 0.99, R1 = 1e-2)
#> Noise DLM block.
#> latent states: 
#>     Noise: Var (1 variable(s))
#> 
#> Linear predictors: 
#>     mu
#> 
#> Status: defined
#> Serie length: 1
#> Interventions at: 
#> Number of latent states: 1
#> Number of linear predictors: 1
```
