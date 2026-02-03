# update_Poisson

Calculate posterior parameter for the Gamma, assuming that the observed
values came from a Poisson model from which the rate parameter (lambda)
have prior distribution Gamma.

## Usage

``` r
update_Poisson(conj.param, ft, Qt, y, parms)
```

## Arguments

- conj.param:

  list: A vector containing the parameters of the Gamma (alpha,beta).

- ft:

  numeric: A vector representing the means from the normal distribution.
  Not used in the default method.

- Qt:

  matrix: A matrix representing the covariance matrix of the normal
  distribution. Not used in the default method.

- y:

  numeric: A vector containing the observations.

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this kernel.

## Value

The parameters of the posterior distribution.

## See also

Other auxiliary functions for a Poisson outcome:
[`convert_Normal_Poisson()`](convert_Normal_Poisson.md),
[`convert_Poisson_Normal()`](convert_Poisson_Normal.md),
[`poisson_pred()`](poisson_pred.md)
