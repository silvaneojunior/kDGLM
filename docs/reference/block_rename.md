# block_rename

block_rename

## Usage

``` r
block_rename(block, pred.names)
```

## Arguments

- block:

  A dlm_block object.

- pred.names:

  A vector of string with names for each linear predictor in block.

## Value

A dlm_block with the linear predictors renamed to the values passed in
names.

## See also

Other auxiliary functions for structural blocks:
[`TF_block()`](tf_block.md), [`block_mult()`](block_mult.md),
[`block_superpos()`](block_superpos.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`intervention()`](intervention.md), [`noise_block()`](noise_block.md),
[`polynomial_block()`](polynomial_block.md),
[`regression_block()`](regression_block.md),
[`specify.dlm_block()`](specify.dlm_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
base.block <- polynomial_block(
  eta = 1,
  order = 1,
  name = "Poly",
  D = 0.95
)

final.block <- block_rename(2 * base.block, c("mu", "sigma"))
```
