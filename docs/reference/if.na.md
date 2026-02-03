# if.na

This function is wrapper for ifelse(is.na(vec),vec,val)

## Usage

``` r
if.na(vec, val)
```

## Arguments

- vec:

  A vector or matrix.

- val:

  The value to replace NA with.

## Value

A vector or matrix with the same dimensions as the input, where any NA
values have been replaced by the specified val argument.
