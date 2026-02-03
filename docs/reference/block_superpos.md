# Auxiliary function for block superposition

An auxiliary function for block superposition.

## Usage

``` r
block_superpos(...)
```

## Arguments

- ...:

  dlm_block: A sequence of block to be combine.

## Value

The combined blocks as a dlm_block.

## Details

Additional details can be found in West and Harrison (1997) , section
6.2.

## References

Mike West, Jeff Harrison (1997). *Bayesian Forecasting and Dynamic
Models (Springer Series in Statistics)*. Springer-Verlag. ISBN
0387947256.

## See also

Other auxiliary functions for structural blocks:
[`TF_block()`](tf_block.md), [`block_mult()`](block_mult.md),
[`block_rename()`](block_rename.md), [`ffs_block()`](ffs_block.md),
[`harmonic_block()`](harmonic_block.md),
[`intervention()`](intervention.md), [`noise_block()`](noise_block.md),
[`polynomial_block()`](polynomial_block.md),
[`regression_block()`](regression_block.md),
[`specify.dlm_block()`](specify.dlm_block.md),
[`summary.dlm_block()`](summary.dlm_block.md)

## Examples

``` r
# Long way
level.1 <- polynomial_block(alpha1 = 1, order = 1)
level.2 <- polynomial_block(alpha2 = 1, order = 2)
season.2 <- harmonic_block(alpha2 = 1, period = 20)

final.block <- block_superpos(level.1, level.2, season.2)

# Short way
final.block <- polynomial_block(alpha1 = 1, order = 1) +
  polynomial_block(alpha2 = 1, order = 2) +
  harmonic_block(alpha2 = 1, period = 20)
```
