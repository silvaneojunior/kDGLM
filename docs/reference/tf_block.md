# Structural blocks for auto regressive trends and regressions

Creates the structure for a Auto Regressive (AR) block (see West and
Harrison (1997) , chapter 9) with desired order. As the package suppose
that the structure of the model is linear, a linearization is applied to
the evolution equation, as described in West and Harrison (1997) ,
chapter 13. This block also supports Transfer Functions, being necessary
to specify the associated pulse when calling the TF_block function (see
arg.).

## Usage

``` r
TF_block(
  ...,
  order,
  noise.var = NULL,
  noise.disc = NULL,
  pulse = 0,
  name = "Var.AR",
  AR.support = "free",
  h = 0,
  a1 = 0,
  R1 = 4,
  monitoring = TRUE,
  multi.states = FALSE,
  D.coef = 1,
  h.coef = 0,
  H.coef = 0,
  a1.coef = c(1, rep(0, order - 1)),
  R1.coef = c(1, rep(0.25, order - 1)),
  monitoring.coef = rep(FALSE, order),
  D.pulse = 1,
  h.pulse = 0,
  H.pulse = 0,
  a1.pulse = 0,
  R1.pulse = 4,
  monitoring.pulse = FALSE
)

AR(
  order = 1,
  noise.var = NULL,
  noise.disc = NULL,
  a1 = 0,
  R1 = 9,
  monitoring = TRUE,
  a1.coef = c(1, rep(0, order - 1)),
  R1.coef = c(1, rep(0.25, order - 1)),
  monitoring.coef = rep(FALSE, order),
  name = "Var.AR",
  X = 1
)

TF(
  pulse,
  order = 1,
  noise.var = NULL,
  noise.disc = NULL,
  a1 = 0,
  R1 = 9,
  monitoring = TRUE,
  a1.coef = c(1, rep(0, order - 1)),
  R1.coef = c(1, rep(0.25, order - 1)),
  monitoring.coef = rep(FALSE, order),
  a1.pulse = 0,
  R1.pulse = 4,
  monitoring.pulse = FALSE,
  name = "Var.AR",
  X = 1
)
```

## Arguments

- ...:

  Named values for the planning matrix.

- order:

  Positive integer: The order of the AR block.

- noise.var:

  Non-negative scalar: The variance of the white noise added to the
  latent state.

- noise.disc:

  Vector or scalar: The value for the discount factor associated with
  the current latent state. If noise.disc is a vector, it should have
  size t and it is interpreted as the discount factor at each observed
  time. If D is a scalar, the same discount will be used for all
  observation.

- pulse:

  Vector or scalar: An optional argument providing the values for the
  pulse for a Transfer Function. Default is 0 (no Transfer Function).

- name:

  String: An optional argument providing the name for this block. Can be
  useful to identify the models with meaningful labels, also, the name
  used will be used in some auxiliary functions.

- AR.support:

  String: Either "constrained" or "free" (default). If AR.support is
  "constrained", then the AR coefficients will be forced to be on the
  interval (-1,1), otherwise, the coefficients will be unrestricted.
  Beware that, under no restriction on the coefficients, there is no
  guarantee that the estimated coefficients will imply in a stationary
  process, furthermore, if the order of the AR block is greater than 1.
  As such the restriction of the coefficients support is only available
  for AR blocks with order equal to 1.

- h:

  Vector or scalar: A drift to be add in the states after the temporal
  evolution (can be interpreted as the mean of the random noise at each
  time). If a vector, it should have size t, and each value will be
  applied in their respective time. If a scalar, the passed value will
  be used for all observations.

- a1:

  Vector or scalar: The prior mean for the states associated with this
  block at time 1. If a1 is a vector, its dimension should be equal to
  the order of the AR block. If a1 is a scalar, its value will be used
  for all coefficients.

- R1:

  Matrix, vector or scalar: The prior covariance matrix for the states
  associated with this block at time 1. If R1 is a matrix, its
  dimensions should be n x n, where n is the order of the AR block. If
  R1 is a vector or scalar, a covariance matrix will be created as a
  diagonal matrix with the values of R1 in the diagonal.

- monitoring:

  bool: A flag indicating if the latent state should be monitored (if
  automated monitoring is used). The default is TRUE.

- multi.states:

  bool: If FALSE (default) a single latent state will be created
  affecting all linear predictor and being affected by all pulses. If
  TRUE, each linear predictor will have its own latent state, but all
  latent states will share the same AR coefficients and all pulse
  effects (each state will have its own pulse though).

- D.coef:

  Array, Matrix, vector or scalar: The values for the discount factors
  associated with the AR coefficients at each time. If D.coef is an
  array, its dimensions should be n x n x t, where n is the order of the
  AR block and t is the length of the outcomes. If D.coef is a matrix,
  its dimensions should be n x n and the same discount matrix will be
  used in all observations. If D.coef is a vector, it should have size t
  and it is interpreted as the discount factor at each observed time
  (same discount for all variable). If D.coef is a scalar, the same
  discount will be used for all AR coefficients at all times.

- h.coef:

  Matrix, vector or scalar: A drift to be add in the AR coefficients
  after the temporal evolution (can be interpreted as the mean of the
  random noise at each time). If a matrix, its dimension should be n x
  t, where n is the order of the AR block and t is the length of the
  series. If a scalar, the passed value will be used for all
  coefficients at each time.

- H.coef:

  Array, Matrix, vector or scalar: The values for the covariance matrix
  for the noise factor associated with the AR coefficients at each time.
  If H.coef is a array, its dimensions should be n x n x t, where n is
  the order of the AR block and t is the length of the outcomes. If
  H.coef is a matrix, its dimensions should be n x n and its values will
  be used for each time. If H.coef is a vector or scalar, a discount
  factor matrix will be created as a diagonal matrix with the values of
  H.coef in the diagonal.

- a1.coef:

  Vector or scalar: The prior mean for the AR coefficients associated
  with this block at time 1. If a1.coef is a vector, its dimension
  should be equal to the order of the AR block. If a1.coef is a scalar,
  its value will be used for all coefficients. If the coefficients are
  restricted to the interval (-1,1), the a1.coef is interpreted as the
  mean for atanh(rho), where rho is the AR coefficient.

- R1.coef:

  Matrix, vector or scalar: The prior covariance matrix for the
  coefficients associated with this block at time 1. If R1.coef is a
  matrix, its dimensions should be n x n, where n is the order of the AR
  block. If R1.coef is a vector or scalar, a covariance matrix will be
  created as a diagonal matrix with the values of R1.coef in the
  diagonal. If the coefficients are restricted to the interval (-1,1),
  the R1.coef is interpreted as the covariance matrix for atanh(rho),
  where rho is the AR coefficient.

- monitoring.coef:

  Vector: A vector of flags indicating which AR coefficients should be
  monitored (if automated monitoring is used). Its size should be n,
  where n is the order of the AR block. The default is that no
  coefficient should be monitored.

- D.pulse:

  Array, Matrix, vector or scalar: The values for the discount factors
  associated with the pulse coefficients at each time. If D.pulse is an
  array, its dimensions should be n x n x t, where n is the number of
  pulses and t is the length of the outcomes. If D.pulse is a matrix,
  its dimensions should be n x n and the same discount matrix will be
  used in all observations. If D.pulse is a vector, it should have size
  t and it is interpreted as the discount factor at each observed time
  (same discount for all variable). If D is a scalar, the same discount
  will be used for all pulse coefficients at all times.

- h.pulse:

  Matrix, vector or scalar: A drift to be add in the pulse effect after
  the temporal evolution (can be interpreted as the mean of the random
  noise at each time). If a matrix, its dimension should be n x t, where
  n is the number of pulses and t is the length of the series. If a
  scalar, the passed value will be used for all latent state at each
  time.

- H.pulse:

  Array, Matrix, vector or scalar: The values for the covariance matrix
  for the noise factor associated with pulse coefficients at each time.
  If H.pulse is an array, its dimensions should be n x n x t, where n is
  the number of pulses and t is the length of the outcomes. If H.pulse
  is a matrix, its dimensions should be n x n and its values will be
  used for each time. If H.pulse is a vector or scalar, a covariance
  matrix will be created as a diagonal matrix with the values of H.pulse
  in the diagonal.

- a1.pulse:

  Vector or scalar: The prior mean for the coefficients associated with
  the pulses at time 1. If a1.pulse is a vector, its dimension should be
  equal to the number of pulses. If a1.pulse is a scalar, its value will
  be used for all coefficients.

- R1.pulse:

  Matrix, vector or scalar: The prior covariance matrix for the
  coefficients associated with the pulses at time 1. If R1.pulse is a
  matrix, its dimensions should be n x n, where n is the number of
  pulses. If R1.pulse is a vector or scalar, a covariance matrix will be
  created as a diagonal matrix with the values of R1.pulse in the
  diagonal.

- monitoring.pulse:

  Vector: A vector of flags indicating which pulse coefficients should
  be monitored (if automated monitoring is used). Its size should be n,
  where n is the number of pulses. The default is that no pulse
  coefficient should be monitored.

- X:

  Vector or scalar: An argument providing the values for the pulse for a
  Transfer Function.

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

- type Character: The type of block (AR).

- AR.support Character: Same as argument.

## Details

For the ..., noise.var, noise.disc, D, H, a1, R1, a1, R1, a1.pulse,
R1.pulse, D.pulse, h.pulse, H.pulse arguments, the user may set one or
more of its values as a string. By doing so, the user will leave the
block partially undefined. The user must then pass the undefined
parameter values as named arguments to the [`fit_model`](fit_model.md)
function. Also, multiple values can be passed, allowing for a
sensitivity analysis for the value of this parameter.

For the details about the implementation see dos Santos et al. (2024) .

For the details about Auto regressive models in the context of DLM's,
see West and Harrison (1997) , chapter 9.

For the details about the linearization of non-linear evolution
equations in the context of DLM's, see West and Harrison (1997) ,
chapter 13.

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
[`block_mult()`](block_mult.md), [`block_rename()`](block_rename.md),
[`block_superpos()`](block_superpos.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`intervention()`](intervention.md), [`noise_block()`](noise_block.md),
[`polynomial_block()`](polynomial_block.md),
[`regression_block()`](regression_block.md),
[`specify.dlm_block()`](specify.dlm_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
#### AR block ####
TF_block(mu = 1, order = 2, noise.disc = 0.9)
#> AR DLM block.
#> latent states: 
#>     Var.AR.State: State.Lag.0, State.Lag.1 (2 variable(s))
#>     Var.AR.Coef: Lag.0, Lag.1 (2 variable(s))
#> 
#> Linear predictors: 
#>     mu
#> 
#> Status: defined
#> Serie length: 1
#> Interventions at: 
#> Number of latent states: 4
#> Number of linear predictors: 1

#### Transfer function ####
TF_block(mu = 1, pulse = beaver1$activ, order = 1, noise.disc = 0.9)
#> Mixed DLM block.
#> latent states: 
#>     Var.AR.State: State.Lag.0 (1 variable(s))
#>     Var.AR.Coef: Lag.0 (1 variable(s))
#>     Var.AR.Pulse.effect: 1 (1 variable(s))
#> 
#> Linear predictors: 
#>     mu
#> 
#> Status: defined
#> Serie length: 114
#> Interventions at: 
#> Number of latent states: 3
#> Number of linear predictors: 1
```
