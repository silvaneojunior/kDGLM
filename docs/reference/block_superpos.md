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
[`TF_block()`](https://silvaneojunior.github.io/kDGLM/reference/tf_block.md),
[`block_mult()`](https://silvaneojunior.github.io/kDGLM/reference/block_mult.md),
[`block_rename()`](https://silvaneojunior.github.io/kDGLM/reference/block_rename.md),
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
level.1 <- polynomial_block(alpha1 = 1, order = 1)
level.2 <- polynomial_block(alpha2 = 1, order = 2)
season.2 <- harmonic_block(alpha2 = 1, period = 20)

final.block <- block_superpos(level.1, level.2, season.2)

# Short way
final.block <- polynomial_block(alpha1 = 1, order = 1) +
  polynomial_block(alpha2 = 1, order = 2) +
  harmonic_block(alpha2 = 1, period = 20)
```
