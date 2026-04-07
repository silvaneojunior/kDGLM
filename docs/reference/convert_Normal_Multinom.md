# convert_Normal_Multinom

Calculate the parameters of the log-Normal that best approximates the
given Dirichlet distribution. The approximation is the best in the sense
that it minimizes the KL divergence from the Dirichlet to the log-Normal

## Usage

``` r
convert_Normal_Multinom(conj.param, parms = list())
```

## Arguments

- conj.param:

  list: A vector containing the concentration parameters of the
  Dirichlet.

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this kernel.

## Value

The parameters of the Normal distribution of the linear predictor.

## See also

Other auxiliary functions for a Multinomial outcome:
[`convert_Multinom_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Multinom_Normal.md),
[`multnom_pred()`](https://silvaneojunior.github.io/kDGLM/reference/multnom_pred.md),
[`update_Multinom()`](https://silvaneojunior.github.io/kDGLM/reference/update_Multinom.md)
