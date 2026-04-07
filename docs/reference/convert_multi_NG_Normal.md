# convert_multi_NG_Normal

Calculate the parameters of the Normal-Gamma that best approximates the
given Multivariate Normal distribution. The distribution obtained for
each outcome is marginal. The approximation is the best in the sense
that it minimizes the KL divergence from the Normal to the Normal-Gamma.
In this approach, we suppose that the first entry of the multivariate
normal represents the mean of the observed data and the second represent
the log variance.

## Usage

``` r
convert_multi_NG_Normal(ft, Qt, parms)
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

Other auxiliary functions for a Normal outcome:
[`multi_normal_gamma_pred()`](https://silvaneojunior.github.io/kDGLM/reference/multi_normal_gamma_pred.md),
[`normal_pred()`](https://silvaneojunior.github.io/kDGLM/reference/normal_pred.md),
[`update_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/update_Normal.md)
