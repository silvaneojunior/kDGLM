# generic_smoother

Generic smoother for all models.

## Usage

``` r
generic_smoother(mt, Ct, at, Rt, G, G.labs, G.idx)
```

## Arguments

- mt:

  matrix: A matrix containing the filtered mean of the latent states at
  each time. Each row should represent one variable.

- Ct:

  array: A 3D-array representing the filtered covariance matrix of the
  latent states at each time. The third dimension should represent the
  time index.

- at:

  matrix: A matrix containing the one-step-ahead mean of the latent
  states at each time based upon the filtered mean. Each row should
  represent one variable.

- Rt:

  array: A 3D-array representing the one-step-ahead covariance matrix of
  the latent states at each time based upon the filtered covariance
  matrix. The third dimension should represent the time index.

- G:

  array: A 3D-array representing the transition matrix of the model at
  each time.

- G.labs:

  matrix: A character matrix containing the type associated with each
  value in G.

- G.idx:

  matrix: A numeric matrix containing the index associated with each
  value in G.

## Value

A list containing the smoothed mean (mts) and covariance (Cts) of the
latent states at each time. Their dimension follows, respectively, the
dimensions of mt and Ct.

## Details

For the models covered in this package, we always assume that the latent
states have Multivariate Normal distribution. With that assumption, we
can use Kalman Smoother algorithm to calculate the posterior of the
states at each time, given everything that has been observed (assuming
that we already know the filtered distribution of the states).

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

[`fit_model`](https://silvaneojunior.github.io/kDGLM/reference/fit_model.md)

[`analytic_filter`](https://silvaneojunior.github.io/kDGLM/reference/analytic_filter.md)
