# convert_Normal_Poisson

Calculate the parameters of the log-Normal that best approximates the
given Gamma distribution. The approximation is the best in the sense
that it minimizes the KL divergence from the Gamma to the log-Normal

## Usage

``` r
convert_Normal_Poisson(conj.param, parms)
```

## Arguments

- conj.param:

  list: A vector containing the parameters of the Gamma (alpha,beta).

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this kernel.

## Value

The parameters of the Normal distribution of the linear predictor.

## See also

Other auxiliary functions for a Poisson outcome:
[`convert_Poisson_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Poisson_Normal.md),
[`poisson_pred()`](https://silvaneojunior.github.io/kDGLM/reference/poisson_pred.md),
[`update_Poisson()`](https://silvaneojunior.github.io/kDGLM/reference/update_Poisson.md)
