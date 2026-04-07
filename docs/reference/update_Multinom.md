# update_Multinom

Calculate posterior parameter for the Dirichlet, assuming that the
observed values came from a Multinomial model from which the number of
trials is known and the prior distribution for the probabilities of each
category have joint distribution Dirichlet.

## Usage

``` r
update_Multinom(conj.param, ft, Qt, y, parms = list())
```

## Arguments

- conj.param:

  list: A vector containing the concentration parameters of the
  Dirichlet.

- ft:

  vector: A vector representing the means from the normal distribution.
  Not used in the default method.

- Qt:

  matrix: A matrix representing the covariance matrix of the normal
  distribution. Not used in the default method.

- y:

  vector: A vector containing the observations.

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this kernel.

## Value

The parameters of the posterior distribution.

## See also

Other auxiliary functions for a Multinomial outcome:
[`convert_Multinom_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Multinom_Normal.md),
[`convert_Normal_Multinom()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Normal_Multinom.md),
[`multnom_pred()`](https://silvaneojunior.github.io/kDGLM/reference/multnom_pred.md)
