# if.null

This function is wrapper for ifelse(is.null(.),.,.)

## Usage

``` r
if.null(vec, val)
```

## Arguments

- vec:

  A vector or matrix.

- val:

  The value to replace NULL with.

## Value

A vector or matrix with the same dimensions as the input, where any NULL
values have been replaced by the specified val argument.
