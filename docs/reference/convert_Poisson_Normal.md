# convert_Poisson_Normal

Calculate the parameters of the Gamma that best approximates the given
log-Normal distribution. The approximation is the best in the sense that
it minimizes the KL divergence from the log-Normal to the Gamma

## Usage

``` r
convert_Poisson_Normal(ft, Qt, parms)
```

## Arguments

- ft:

  numeric: A vector representing the means from the normal distribution.

- Qt:

  matrix: A matrix representing the covariance matrix of the normal
  distribution.

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this kernel.

## Value

The parameters of the conjugated distribution of the linear predictor.

## See also

Other auxiliary functions for a Poisson outcome:
[`convert_Normal_Poisson()`](convert_Normal_Poisson.md),
[`poisson_pred()`](poisson_pred.md),
[`update_Poisson()`](update_Poisson.md)
