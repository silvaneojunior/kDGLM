# array_mult_left

Calculates the matrix product between an array and a matrix.

## Usage

``` r
array_mult_left(A, B)
```

## Arguments

- A:

  A 3-D array with shapes n x m x k.

- B:

  A matrix with shapes m x l.

## Details

For an array A with shapes n x m x k and a matrix B with shape m x l,
this operations returns an array C, with shapes n x l x k, so that
C\[,,i\] = A\[,,i\]
