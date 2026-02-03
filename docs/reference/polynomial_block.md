# Structural blocks for polynomial trends and regressions

Creates the structure for a polynomial block with desired order.

## Usage

``` r
polynomial_block(
  ...,
  order = 1,
  name = "Var.Poly",
  D = 1,
  h = 0,
  H = 0,
  a1 = 0,
  R1 = c(9, rep(1, order - 1)),
  monitoring = c(TRUE, rep(FALSE, order - 1))
)

pol(
  order = 1,
  D = 0.95,
  a1 = 0,
  R1 = 9,
  monitoring = c(TRUE, rep(FALSE, order - 1)),
  name = "Var.Poly",
  X = 1
)
```

## Arguments

- ...:

  Named values for the planning matrix.

- order:

  Positive integer: The order of the polynomial structure.

- name:

  String: An optional argument providing the name for this block. Can be
  useful to identify the models with meaningful labels, also, the name
  used will be used in some auxiliary functions.

- D:

  Array, Matrix, vector or scalar: The values for the discount factors
  associated with the latent states at each time. If D is an array, its
  dimensions should be n x n x t, where n is the order of the polynomial
  block and t is the length of the outcomes. If D is a matrix, its
  dimensions should be n x n and the same discount matrix will be used
  in all observations. If D is a vector, it should have size t and it is
  interpreted as the discount factor at each observed time (same
  discount for all variable). If D is a scalar, the same discount will
  be used for all latent states at all times.

- h:

  Matrix, vector or scalar: A drift to be add after the temporal
  evolution (can be interpreted as the mean of the random noise at each
  time). If a matrix, its dimension should be n x t, where n is the
  number of latent states (i.e., the order) and t is the length of the
  series. If a vector, it should have size t, and each value will be
  applied to the first latent state (the one which affects the linear
  predictors) in their respective time. If a scalar, the passed value
  will be used for the first latent state at each time.

- H:

  Array, Matrix, vector or scalar: The values for the covariance matrix
  for the noise factor at each time. If H is an array, its dimensions
  should be n x n x t, where n is the order of the polynomial block and
  t is the length of the series. If H is a matrix, its dimensions should
  be n x n and its values will be used for each time. If H is a vector
  or scalar, a discount factor matrix will be created as a diagonal
  matrix with the values of H in the diagonal.

- a1:

  Vector or scalar: The prior mean for the latent states associated with
  this block at time 1. If a1 is a vector, its dimension should be equal
  to the order of the polynomial block. If a1 is a scalar, its value
  will be used for all latent states.

- R1:

  Matrix, vector or scalar: The prior covariance matrix for the latent
  states associated with this block at time 1. If R1 is a matrix, its
  dimensions should be n x n. If R1 is a vector or scalar, a covariance
  matrix will be created as a diagonal matrix with the values of R1 in
  the diagonal.

- monitoring:

  Vector: A vector of flags indicating which variables should be
  monitored (if automated monitoring is used). Its size should be n. The
  default is that only the first order component of this structure
  should be monitored.

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

- G.idx Matrix: A n x n character matrix containing the index each
  element of G.

- D Array: A 3D-array containing the discount factor matrix for each
  time. Its dimension should be n x n x t, where n is the number of
  latent states and t is the time series length.

- h Matrix: The mean for the random noise of the temporal evolution. Its
  dimension should be n x t.

- H Array: A 3D-array containing the covariance matrix of the noise for
  each time. Its dimension should be the same as D.

- a1 Vector: The prior mean for the latent vector.

- R1 Matrix: The prior covariance matrix for the latent vector.

- var.names list: A list containing the variables indexes by their name.

- order Positive integer: Same as argument.

- n Positive integer: The number of latent states associated with this
  block (same value as order).

- t Positive integer: The number of time steps associated with this
  block. If 1, the block is compatible with blocks of any time length,
  but if t is greater than 1, this block can only be used with blocks of
  the same time length.

- k Positive integer: The number of outcomes associated with this block.
  This block can only be used with blocks with the same outcome length.

- pred.names Vector: The name of the linear predictors associated with
  this block.

- monitoring Vector: Same as argument.

- type Character: The type of block (polynomial).

## Details

For the ..., D, H, a1 and R1 arguments, the user may set one or more of
its values as a string. By doing so, the user will leave the block
partially undefined. The user must then pass the undefined parameter
values as named arguments to the [`fit_model`](fit_model.md) function.
Also, multiple values can be passed, allowing for a sensitivity analysis
for the value of this parameter.

For the details about the implementation see dos Santos et al. (2024) .

For the details about polynomial trend in the context of DLM's, see West
and Harrison (1997) , chapter 7.

For the details about dynamic regression models in the context of DLM's,
see West and Harrison (1997) , chapters 6 and 9.

## References

Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”  
  
Mike West, Jeff Harrison (1997). *Bayesian Forecasting and Dynamic
Models (Springer Series in Statistics)*. Springer-Verlag. ISBN
0387947256.

## See also

[`fit_model`](fit_model.md)

Other auxiliary functions for structural blocks:
[`TF_block()`](tf_block.md), [`block_mult()`](block_mult.md),
[`block_rename()`](block_rename.md),
[`block_superpos()`](block_superpos.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`intervention()`](intervention.md), [`noise_block()`](noise_block.md),
[`regression_block()`](regression_block.md),
[`specify.dlm_block()`](specify.dlm_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
# Creating a first order structure for a model with 2 outcomes.
# One block is created for each outcome
# with each block being associated with only one of the outcomes.
level.1 <- polynomial_block(alpha1 = 1, order = 1)
level.2 <- polynomial_block(alpha2 = 1, order = 1)

# Creating a block with shared effect between the outcomes
level.3 <- polynomial_block(alpha1 = 1, alpha2 = 1, order = 2)
```
