# update_NG

Calculate posterior parameter for the Normal-Gamma, assuming that the
observed values came from a Normal model from which the prior
distribution for the mean and the precision have joint distribution
Normal-Gamma

## Usage

``` r
update_NG(conj.param, ft, Qt, y, parms = list())
```

## Arguments

- conj.param:

  list: A vector containing the parameters of the Normal-Gamma
  (mu0,c0,alpha,beta).

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
