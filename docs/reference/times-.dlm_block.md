# \*.fitted_dlm

Define product operator for class dlm_block. This method is wrapper for
the block_mult function.

## Usage

``` r
# S3 method for class 'dlm_block'
e1 * e2
```

## Arguments

- e1:

  A dlm_block (if e2 is an integer) or an integer (if e2 is a
  dlm_block).

- e2:

  An integer (if e1 is an dlm_block) or a dlm_block (if e1 is an
  integer).

## Value

The combined replicated blocks as a dlm_block.

## See also

[`block_mult`](block_mult.md)
