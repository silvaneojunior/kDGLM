# update_Gamma

Calculate posterior parameter for the Inverse-Gamma, assuming that the
observed values came from a Gamma model from which the shape parameter
(phi) is known and the mean (mu) have prior distribution Inverse-Gamma.

## Usage

``` r
update_Gamma(conj.param, ft, Qt, y, parms)
```

## Arguments

- conj.param:

  list: A vector containing the parameters of the Inverse-Gamma
  (alpha,beta).

- ft:

  numeric: A vector representing the means from the normal distribution.
  Not used in the default method.

- Qt:

  matrix: A matrix representing the covariance matrix of the normal
  distribution. Not used in the default method.

- y:

  numeric: A vector containing the observations.

- parms:

  list: A list of extra known parameters of the distribution. For this
  kernel, parms should containing the shape parameter (phi) for the
  observational gamma model.

## Value

The parameters of the posterior distribution.

## See also

Other auxiliary functions for a Gamma outcome with known shape:
[`convert_Gamma_Normal()`](convert_Gamma_Normal.md),
[`convert_Normal_Gamma()`](convert_Normal_Gamma.md),
[`gamma_pred()`](gamma_pred.md)
