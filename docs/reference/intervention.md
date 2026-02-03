# An auxiliary function for model intervention

This function adds timely modifications to a dlm_block, such that in the
specified time the model will override the usual value of the each
variable to the value chosen by the user.

## Usage

``` r
intervention(
  block,
  time,
  var.index = 1:block$n,
  FF = NULL,
  D = NULL,
  h = NULL,
  H = NULL,
  G = NULL
)
```

## Arguments

- block:

  dlm_block: The block to add the intervention.

- time:

  Vector: A sequence of integers indicating the time of the
  intervention.

- var.index:

  Vector: A sequence of integers indicating which variables should be
  modified in the intervention.

- FF:

  Array: A n x k x t array with the modified FF to be used during the
  intervention, where n is the length of var.index, k is the number of
  linear predictors in the block and t is the size of time (can be
  omitted if time is a scalar).

- D:

  Array: A n x n x t array with the modified D to be used during the
  intervention, where n is the length of var.index and t is the size of
  time (can be omitted if time is a scalar).

- h:

  matrix: A n x t matrix with the modified h to be used during the
  intervention, where n is the length of var.index and t is the size of
  time (can be omitted if time is a scalar).

- H:

  Array: A n x n x t array with the modified H to be used during the
  intervention, where n is the length of var.index and t is the size of
  time (can be omitted if time is a scalar).

- G:

  Array: A n x n x t array with the modified G to be used during the
  intervention, where n is the length of var.index and t is the size of
  time (can be omitted if time is a scalar).

## Value

A dlm_block with the added intervention.

## See also

Other auxiliary functions for structural blocks:
[`TF_block()`](tf_block.md), [`block_mult()`](block_mult.md),
[`block_rename()`](block_rename.md),
[`block_superpos()`](block_superpos.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`noise_block()`](noise_block.md),
[`polynomial_block()`](polynomial_block.md),
[`regression_block()`](regression_block.md),
[`specify.dlm_block()`](specify.dlm_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
data <- c(AirPassengers)
# Adding an artificial change, so that we can make an intervention on the data at that point
# Obviously, one should NOT change their own data.
data[60:144] <- data[60:144] + 500

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

# Reducing the discount factor so that the model can capture the expected change.
level <- level |> intervention(time = 60, H = 1, var.index = 1)
# Comment the line above to see the fit without the intervention

outcome <- Poisson(lambda = "rate", data = data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)

plot(fitted.data, plot.pkg = "base")

```
