# rmvnorm

Obtains a sample from a multivariate normal distribution.

## Usage

``` r
rmvnorm(n, mu, Sigma, norm.x = matrnorm(k, n, seed = round(runif(1) * 1e+15)))
```

## Arguments

- n:

  integer: The sample size.

- mu:

  numeric: The mean vector

- Sigma:

  matrix: The Covariance matrix.
