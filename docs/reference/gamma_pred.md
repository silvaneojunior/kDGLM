# gamma_pred

Calculate the values for the predictive distribution given the values of
the parameter of the conjugated distribution of the linear predictor.
The data is assumed to have Gamma distribution with known shape
parameter phi and its mean having distribution Inverse-Gamma with shape
parameter a e rate parameter b. In this scenario, the marginal
distribution of the data is Beta prime with parameters phi, alpha, beta
/ phi.

## Usage

``` r
gamma_pred(conj.param, outcome = NULL, parms = list(), pred.cred = 0.95)
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
  must contain the shape parameter phi of the observational model.

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

Other auxiliary functions for a Gamma outcome with known shape:
[`convert_Gamma_Normal()`](convert_Gamma_Normal.md),
[`convert_Normal_Gamma()`](convert_Normal_Gamma.md),
[`update_Gamma()`](update_Gamma.md)
