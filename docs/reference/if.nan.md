# if.nan

This function is wrapper for ifelse(is.nan(vec),vec,val)

## Usage

``` r
if.nan(vec, val)
```

## Arguments

- vec:

  A vector or matrix.

- val:

  The value to replace NaN with.

## Value

A vector or matrix with the same dimensions as the input, where any NaN
values have been replaced by the specified val argument.
