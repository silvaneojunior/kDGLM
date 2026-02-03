# Specify method for dlm blocks

Sets the values of undefined parameters in a block to those passed by
the user.

## Usage

``` r
# S3 method for class 'dlm_block'
specify(x, ...)
```

## Arguments

- x:

  dlm_block: A undefined dlm_block object from which the undefined
  parameters shall be substituted.

- ...:

  A set of named values for each unknown parameter.

## Value

The initual block, but with the undefined parameters set to the chosen
values.

## See also

Other auxiliary functions for structural blocks:
[`TF_block()`](tf_block.md), [`block_mult()`](block_mult.md),
[`block_rename()`](block_rename.md),
[`block_superpos()`](block_superpos.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`intervention()`](intervention.md), [`noise_block()`](noise_block.md),
[`polynomial_block()`](polynomial_block.md),
[`regression_block()`](regression_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
season <- harmonic_block(rate = 1, period = 12, D = "D.sazo") |>
  specify(D.sazo = 0.975)
```
