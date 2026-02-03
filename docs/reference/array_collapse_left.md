# array_collapse_left

Calculates the matrix product between an array and a vector.

## Usage

``` r
array_collapse_left(A, B)
```

## Arguments

- A:

  A 3-D array with shapes n x m x k.

- B:

  A matrix with shapes m x 1.

## Details

For an array A with shapes n x m x k and a vector B with shape m, this
operations returns a matrix C, with shapes n x k, so that C\[,i\] =
A\[,,i\]
