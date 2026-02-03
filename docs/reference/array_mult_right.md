# array_mult_right

Calculates the matrix product between an array and a matrix.

## Usage

``` r
array_mult_right(A, B)
```

## Arguments

- A:

  A 3-D array with shapes n x m x k.

- B:

  A matrix with shapes l x n.

## Details

For an array A with shapes m x n x k and a matrix B with shape l x m,
this operations returns an array C, with shapes l x n x k, so that
C\[,,i\] = B
