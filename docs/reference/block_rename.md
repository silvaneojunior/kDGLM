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
[`TF_block()`](https://silvaneojunior.github.io/kDGLM/reference/tf_block.md),
[`block_mult()`](https://silvaneojunior.github.io/kDGLM/reference/block_mult.md),
[`block_superpos()`](https://silvaneojunior.github.io/kDGLM/reference/block_superpos.md),
[`ffs_block()`](https://silvaneojunior.github.io/kDGLM/reference/ffs_block.md),
[`harmonic_block()`](https://silvaneojunior.github.io/kDGLM/reference/harmonic_block.md),
[`intervention()`](https://silvaneojunior.github.io/kDGLM/reference/intervention.md),
[`noise_block()`](https://silvaneojunior.github.io/kDGLM/reference/noise_block.md),
[`polynomial_block()`](https://silvaneojunior.github.io/kDGLM/reference/polynomial_block.md),
[`regression_block()`](https://silvaneojunior.github.io/kDGLM/reference/regression_block.md),
[`specify.dlm_block()`](https://silvaneojunior.github.io/kDGLM/reference/specify.dlm_block.md),
[`summary.dlm_block()`](https://silvaneojunior.github.io/kDGLM/reference/summary.dlm_block.md)

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
