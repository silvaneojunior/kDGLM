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
[`TF_block()`](tf_block.md), [`block_rename()`](block_rename.md),
[`block_superpos()`](block_superpos.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`intervention()`](intervention.md), [`noise_block()`](noise_block.md),
[`polynomial_block()`](polynomial_block.md),
[`regression_block()`](regression_block.md),
[`specify.dlm_block()`](specify.dlm_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
# Long way
level <- polynomial_block(alpha = 1, order = 1)

final.block <- block_mult(level, 5)

# Short way
final.block <- 5 * polynomial_block(alpha = 1, order = 1)
```
