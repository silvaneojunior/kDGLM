# Auxiliary function to replicate blocks

An auxiliary function to replicate blocks.

## Usage

``` r
block_mult(block, k)
```

## Arguments

- block:

  dlm_block: A block to be replicated

- k:

  Integer: The number of blocks to generate.

## Value

The combined replicated blocks as a dlm_block.

## See also

Other auxiliary functions for structural blocks:
[`TF_block()`](https://silvaneojunior.github.io/kDGLM/reference/tf_block.md),
[`block_rename()`](https://silvaneojunior.github.io/kDGLM/reference/block_rename.md),
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
# Long way
level <- polynomial_block(alpha = 1, order = 1)

final.block <- block_mult(level, 5)

# Short way
final.block <- 5 * polynomial_block(alpha = 1, order = 1)
```
