# analytic_filter

Fit a model given the observed value and the model parameters.

## Usage

``` r
analytic_filter(
  outcomes,
  a1 = 0,
  R1 = 1,
  FF,
  FF.labs,
  G,
  G.labs,
  G.idx,
  D,
  h,
  H,
  p.monit = NA,
  monitoring = FALSE,
  safe.mode = TRUE
)
```

## Arguments

- outcomes:

  list: The observed data. It should contain objects of the class
  dlm_distr.

- a1:

  numeric: The prior mean at the latent vector.

- R1:

  matrix: The prior covariance matrix at the latent vector.

- FF:

  array: A 3D-array containing the planning matrix at each time. Its
  dimension should be n x k x t, where n is the number of latent states,
  k is the number of linear predictors in the model and t is the time
  series length.

- FF.labs:

  matrix: A character matrix containing the label associated with each
  value in FF.

- G:

  array: A 3D-array containing the evolution matrix at each time. Its
  dimension should be n x n x t, where n is the number of latent states
  and t is the time series length.

- G.labs:

  matrix: A character matrix containing the label associated with each
  value in G.

- G.idx:

  matrix: A numeric matrix containing the index associated with each
  value in G.

- D:

  array: A 3D-array containing the discount factor matrix at each time.
  Its dimension should be n x n x t, where n is the number of latent
  states and t is the time series length.

- h:

  matrix: A drift to be added after the temporal evolution (can be
  interpreted as the mean of the random noise at each time). Its
  dimension should be n x t, where t is the length of the series and n
  is the number of latent states.

- H:

  array: A 3D-array containing the covariance matrix of the noise at
  each time. Its dimension should be the same as D.

- p.monit:

  numeric (optional): The prior probability of changes in the latent
  space variables that are not part of its dynamic.

- monitoring:

  numeric: A vector of flags indicating which latent states should be
  monitored.

- safe.mode:

  boolean: A flag indicating if consistency check should be performed at
  each time step. Recommended to be left on, but if you know what you
  are doing (i.e., you tested the model and it is safe) and need to fit
  it several times, you can disable the checks to save some time.

## Value

A list containing the following values:

- mt matrix: The filtered mean of the latent states for each time.
  Dimensions are n x t.

- Ct array: A 3D-array containing the filtered covariance matrix of the
  latent states for each time. Dimensions are n x n x t.

- at matrix: The one-step-ahead mean of the latent states at each time.
  Dimensions are n x t.

- Rt array: A 3D-array containing the one-step-ahead covariance matrix
  for latent states at each time. Dimensions are n x n x t.

- ft matrix: The one-step-ahead mean of the linear predictors at each
  time. Dimensions are k x t.

- Qt array: A 3D-array containing the one-step-ahead covariance matrix
  for linear predictors at each time. Dimensions are k x k x t.

- ft.star matrix: The filtered mean of the linear predictors for each
  time. Dimensions are k x t.

- Qt.star array: A 3D-array containing the linear predictors matrix of
  the latent state for each time. Dimensions are k x k x t.

- FF array: The same as the argument (same values).

- G matrix: The same as the argument (same values).

- G.labs matrix: The same as the argument (same values).

- G.idx matrix: The same as the argument (same values).

- D array: The same as the argument (same values).

- h array: The same as the argument (same values).

- H array: The same as the argument (same values).

- W array: A 3D-array containing the effective covariance matrix of the
  noise for each time, i.e., considering both H and D. Its dimension are
  the same as H and D.

- monitoring numeric: The same as the argument (same values).

- outcomes list: The same as the argument outcomes (same values).

- pred.names numeric: The names of the linear predictors.

- safe.mode bool: The same as the argument outcomes (same values).

## Details

For the models covered in this package, we always use the approach
described in Alves et al. (2024) , including, in particular, the
filtering algorithm presented in that work.

For the details about the implementation see dos Santos et al. (2024) .

For the details about the algorithm implemented see Alves et al. (2024)
, Petris et al. (2009) , chapter 2, West and Harrison (1997) , chapter
4, and Kalman (1960) .

## References

Mariane Branco Alves, Helio S. Migon, Raíra Marotta, Junior, Silvaneo
Vieira dos Santos (2024). “k-parametric Dynamic Generalized Linear
Models: a sequential approach via Information Geometry.” 2201.05387.  
  
Junior, Silvaneo Vieira dos Santos, Mariane Branco Alves, Helio S. Migon
(2024). “kDGLM: an R package for Bayesian analysis of Dynamic
Generialized Linear Models.”  
  
Rudolph Emil Kalman (1960). “A New Approach to Linear Filtering and
Prediction Problems.” *Transactions of the ASME–Journal of Basic
Engineering*, **82**(Series D), 35–45.  
  
Giovanni Petris, Sonia Petrone, Patrizia Campagnoli (2009). *Dynamic
Linear Models with R*, useR! Springer-Verlag, New York.  
  
Mike West, Jeff Harrison (1997). *Bayesian Forecasting and Dynamic
Models (Springer Series in Statistics)*. Springer-Verlag. ISBN
0387947256.

## See also

[`fit_model`](fit_model.md)

[`generic_smoother`](generic_smoother.md)
