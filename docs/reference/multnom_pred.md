# multnom_pred

Calculate the values for the predictive distribution given the values of
the parameter of the conjugated distribution of the linear predictor.
The data is assumed to have Multinomial distribution with known number
of trial N and the probability vector having distribution Dirichlet with
parameters alpha_i. In this scenario, the marginal distribution of the
data is Dirichlet-Multinomial with parameters N and alpha_i.

## Usage

``` r
multnom_pred(conj.param, outcome, parms = list(), pred.cred = 0.95)
```

## Arguments

- conj.param:

  List or data.frame: The parameters of the conjugated distributions of
  the linear predictor.

- outcome:

  Vector or matrix: The observed values at the current time. The value
  passed is used to compute N.

- parms:

  List (optional): A list of extra parameters for the model. Not used in
  this function.

- pred.cred:

  Numeric: the desired credibility for the credibility interval.

## Value

A list containing the following values:

- pred vector/matrix: the mean of the predictive distribution of a next
  observation. Same type and shape as the parameter in model.

- var.pred vector/matrix: the variance of the predictive distribution of
  a next observation. Same type and shape as the parameter in model.

- icl.pred vector/matrix: the percentile of 100\*((1-pred.cred)/2)

- icu.pred vector/matrix: the percentile of 100\*(1-(1-pred.cred)/2)

- log.like vector: the The log likelihood for the outcome given the
  conjugated parameters.

## See also

Other auxiliary functions for a Multinomial outcome:
[`convert_Multinom_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Multinom_Normal.md),
[`convert_Normal_Multinom()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Normal_Multinom.md),
[`update_Multinom()`](https://silvaneojunior.github.io/kDGLM/reference/update_Multinom.md)
