# poisson_pred

Calculate the values for the predictive distribution given the values of
the parameter of the conjugated distribution of the linear predictor.
The data is assumed to have Poisson distribution with its mean having
distribution Gamma with shape parameter a e rate parameter b. In this
scenario, the marginal distribution of the data is Negative Binomial
with a as the dispersion parameter and b/(b+1) as the probability.

## Usage

``` r
poisson_pred(conj.param, outcome = NULL, parms = list(), pred.cred = 0.95)
```

## Arguments

- conj.param:

  list or data.frame: The parameters of the conjugated distributions of
  the linear predictor.

- outcome:

  numeric or matrix (optional): The observed values at the current time.
  Not used in this function.

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

Other auxiliary functions for a Poisson outcome:
[`convert_Normal_Poisson()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Normal_Poisson.md),
[`convert_Poisson_Normal()`](https://silvaneojunior.github.io/kDGLM/reference/convert_Poisson_Normal.md),
[`update_Poisson()`](https://silvaneojunior.github.io/kDGLM/reference/update_Poisson.md)
