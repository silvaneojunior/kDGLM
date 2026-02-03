# Basic structural blocks

Creates the basic structure for a dlm block with desired order.

## Usage

``` r
base_block(..., order, name, D, h, H, a1, R1, monitoring)
```

## Arguments

- ...:

  Named values for the planning matrix.

- order:

  integer: The order of the structure. Must be positive

- name:

  character: An optional argument providing the name for this block. Can
  be useful to identify the models with meaningful labels, also, the
  name used will be used in some auxiliary functions.

- D:

  array, matrix, vector or scalar: The values for the discount factors
  associated with the latent states at each time. If D is an array, its
  dimensions should be n x n x t, where n is the order of the polynomial
  block and t is the length of the outcomes. If D is a matrix, its
  dimensions should be n x n and the same discount matrix will be used
  in all observations. If D is a vector, it should have size t and it is
  interpreted as the discount factor at each observed time (same
  discount for all variable). If D is a scalar, the same discount will
  be used for all latent states at all times.

- h:

  matrix, vector or scalar: A drift to be add after the temporal
  evolution (can be interpreted as the mean of the random noise at each
  time). If a matrix, its dimension should be n x t, where n is the
  number of latent states (i.e., the order) and t is the length of the
  series. If a vector, it should have size t, and each value will be
  applied to the first latent state (the one which affects the linear
  predictors) in their respective time. If a scalar, the passed value
  will be used for the first latent state at each time.

- H:

  array, matrix, vector or scalar: The values for the covariance matrix
  for the noise factor at each time. If H is an array, its dimensions
  should be n x n x t, where n is the order of the polynomial block and
  t is the length of the series. If H is a matrix, its dimensions should
  be n x n and its values will be used for each time. If H is a vector
  or scalar, a discount factor matrix will be created as a diagonal
  matrix with the values of H in the diagonal.

- a1:

  vector or scalar: The prior mean for the latent states associated with
  this block at time 1. If a1 is a vector, its dimension should be equal
  to the order of the polynomial block. If a1 is a scalar, its value
  will be used for all latent states.

- R1:

  matrix, vector or scalar: The prior covariance matrix for the latent
  states associated with this block at time 1. If R1 is a matrix, its
  dimensions should be n x n. If R1 is a vector or scalar, a covariance
  matrix will be created as a diagonal matrix with the values of R1 in
  the diagonal.

- monitoring:

  vector: A vector of flags indicating which variables should be
  monitored (if automated monitoring is used). Its size should be n. The
  default is that only the first order component of this structure
  should be monitored.
