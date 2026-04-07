# update_Normal

Calculate posterior parameter for the Normal, assuming that the observed
values came from a Normal model from which the covariance is known and
the prior distribution for the mean vector have Normal distribution

## Usage

``` r
update_Normal(conj.param, ft, Qt, y, parms)
```

## Arguments

- conj.param:

  list: A vector containing the concentration parameters of the Normal.

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
  kernel, parms should containing the covariance matrix parameter (V)
  for the observational Normal model.

## Value

The parameters of the posterior distribution.

## See also

Other auxiliary functions for a Normal outcome:
[`convert_multi_NG_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_multi_NG_Normal.md),
[`multi_normal_gamma_pred()`](https://silvaneojunior.github.io/kDGLM/reference/multi_normal_gamma_pred.md),
[`normal_pred()`](https://silvaneojunior.github.io/kDGLM/reference/normal_pred.md)
