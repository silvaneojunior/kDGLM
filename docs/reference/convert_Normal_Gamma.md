# convert_Normal_Gamma

Calculates the parameters of the log-Normal that best approximates the
given Inverse-Gamma distribution. The approximation is the best in the
sense that it minimizes the KL divergence from the Inverse-Gamma to the
log-Normal

## Usage

``` r
convert_Normal_Gamma(conj.param, parms)
```

## Arguments

- conj.param:

  list: A vector containing the parameters of the Inverse-Gamma
  (alpha,beta).

- parms:

  list: A list of extra known parameters of the distribution. Not used
  in this function.

## Value

The parameters of the Normal distribution of the linear predictor.

## See also

Other auxiliary functions for a Gamma outcome with known shape:
[`convert_Gamma_Normal()`](convert_Gamma_Normal.md),
[`gamma_pred()`](gamma_pred.md), [`update_Gamma()`](update_Gamma.md)
