# normal_pred

Calculate the values for the predictive distribution given the values of
the parameter of the conjugated distribution of the linear predictor.
The data is assumed to have Normal distribution with known variance and
its mean having distribution Normal. In this scenario, the marginal
distribution of the data is also Normal.

## Usage

``` r
normal_pred(conj.param, outcome = NULL, parms = list(), pred.cred = 0.95)
```

## Arguments

- conj.param:

  list or data.frame: The parameters of the conjugated distributions of
  the linear predictor.

- outcome:

  numeric or matrix (optional): The observed values at the current time.
  Not used in this function.

- parms:

  list: A list of extra parameters for the model. For this function, it
  must contain the observational covariance matrix, V

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
[`convert_multi_NG_Normal()`](convert_multi_NG_Normal.md),
[`multi_normal_gamma_pred()`](multi_normal_gamma_pred.md),
[`update_Normal()`](update_Normal.md)
