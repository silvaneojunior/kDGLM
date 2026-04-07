# convert_Gamma_Normal

Calculate the parameters of the Inverse-Gamma that best approximates the
given log-Normal distribution. The approximation is the best in the
sense that it minimizes the KL divergence from the log-Normal to the
Inverse-Gamma

## Usage

``` r
convert_Gamma_Normal(ft, Qt, parms)
```

## Arguments

- ft:

  vector: A vector representing the means from the normal distribution.

- Qt:

  matrix: A matrix representing the covariance matrix of the normal
  distribution.

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this function.

## Value

The parameters of the conjugated distribution of the linear predictor.

## See also

Other auxiliary functions for a Gamma outcome with known shape:
[`convert_Normal_Gamma()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Normal_Gamma.md),
[`gamma_pred()`](https://silvaneojunior.github.io/kDGLM/reference/gamma_pred.md),
[`update_Gamma()`](https://silvaneojunior.github.io/kDGLM/reference/update_Gamma.md)
