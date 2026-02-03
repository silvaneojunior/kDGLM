# convert_Multinom_Normal

Calculate the parameters of the Dirichlet that best approximates the
given log-Normal distribution. The approximation is the best in the
sense that it minimizes the KL divergence from the log-Normal to the
Dirichlet.

## Usage

``` r
convert_Multinom_Normal(ft, Qt, parms = list())
```

## Arguments

- ft:

  vector: A vector representing the means from the normal distribution.

- Qt:

  matrix: A matrix representing the covariance matrix of the normal
  distribution.

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this kernel.

## Value

The parameters of the conjugated distribution of the linear predictor.

## See also

Other auxiliary functions for a Multinomial outcome:
[`convert_Normal_Multinom()`](convert_Normal_Multinom.md),
[`multnom_pred()`](multnom_pred.md),
[`update_Multinom()`](update_Multinom.md)
