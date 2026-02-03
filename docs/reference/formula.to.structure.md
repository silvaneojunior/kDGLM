# formula.to.structure

formula.to.structure

## Usage

``` r
# S3 method for class 'to.structure'
formula(formula, data, label = "mu")
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted.

- data:

  an optional data frame, list or environment (or object coercible by
  as.data.frame to a data frame) containing the variables in the model.
  If not found in data, the variables are taken from
  environment(formula), typically the environment from which glm is
  called.

- label:

  An optional character naming the linear predictor.
