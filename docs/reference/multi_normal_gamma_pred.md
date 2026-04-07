# update_multi_NG_chol multi_normal_gamma_pred

update_multi_NG_chol multi_normal_gamma_pred

## Usage

``` r
multi_normal_gamma_pred(
  conj.param,
  outcome = NULL,
  parms = list(),
  pred.cred = 0.95
)
```

## Arguments

- conj.param:

  list or data.frame: The parameters of the conjugated distribution
  (Normal-Gamma) of the linear predictor.

- outcome:

  numeric or matrix (optional): The observed values at the current time.

- parms:

  list (optional): A list of extra parameters for the model. Not used in
  this function.

- pred.cred:

  numeric: the desired credibility for the credibility interval.

## Value

A list containing the following values:

- pred numeric/matrix: the mean of the predictive distribution of a next
  observation. Same type and shape as the parameter in model.

- var.pred numeric/matrix: the variance of the predictive distribution
  of a next observation. Same type and shape as the parameter in model.

- icl.pred numeric/matrix: the percentile of 100\*((1-pred.cred)/2)

- icu.pred numeric/matrix: the percentile of 100\*(1-(1-pred.cred)/2)

- log.like numeric: the The log likelihood for the outcome given the
  conjugated parameters.

## See also

Other auxiliary functions for a Normal outcome:
[`convert_multi_NG_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_multi_NG_Normal.md),
[`normal_pred()`](https://silvaneojunior.github.io/kDGLM/reference/normal_pred.md),
[`update_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/update_Normal.md)
