# Fitting kDGLM models

Fit a model given its structure and the observed data. This function can
be used for any supported family (see vignette).

## Usage

``` r
fit_model(
  ...,
  smooth = TRUE,
  p.monit = NA,
  condition = "TRUE",
  metric = "log.like",
  lag = 1,
  pred.cred = 0.95,
  metric.cutoff = NA,
  save.models = FALSE,
  silent = FALSE,
  safe.mode = TRUE
)
```

## Arguments

- ...:

  dlm_block or dlm_distr objects or named values: The structural blocks
  of the model (dlm_block objects), alongside the model outcomes
  (dlm_distr object). If at least one block is undefined, the user must
  also provide its value in this argument (see last example).

- smooth:

  boolean: A flag indicating if the smoothed distribution for the latent
  states should be calculated.

- p.monit:

  numeric (optional): The prior probability of changes in the latent
  space variables that are not part of its dynamic. Only used when
  performing sensitivity analysis.

- condition:

  character: A character defining which combinations of undefined hyper
  parameter should be tested. See example for details. Only used when
  performing sensitivity analysis.

- metric:

  character: The name of the metric to use for model selection. One of
  log-likelihood for the one-step-ahead prediction ("log.like"), Mean
  Absolute Scaled Error ("mase") (Hyndman and Koehler 2006) or Interval
  Score ("interval.score") (Bracher et al. 2021) . Only used when
  performing sensitivity analysis.

- lag:

  integer: The number of steps ahead used for the prediction when
  calculating the metrics. If lag\<0, predictions are made using the
  smoothed distribution of the latent states. Only used when performing
  sensitivity analysis.

- pred.cred:

  numeric: A number between 0 and 1 (not included) indicating the
  credibility interval for predictions. If not within the valid range of
  values, 0.95 will be used. Only used when performing sensitivity
  analysis.

- metric.cutoff:

  integer: The number of observations to ignore when calculating the
  metrics. Default is 1/10 of the number of observations (rounded down).
  Only used when performing sensitivity analysis.

- save.models:

  boolean: A flag indicating if all evaluated models should be saved. If
  FALSE, only the best model (according to the chosen metric) will be
  saved. Only used when performing sensitivity analysis.

- silent:

  boolean: A flag indicating if a progress bar should be printed. Only
  used when performing sensitivity analysis.

- safe.mode:

  boolean: A flag indicating if consistency check should be performed at
  each time step. Recommended to be left on, but if you know what you
  are doing (i.e., you tested the model and it is safe) and need to fit
  it several times, you can disable the checks to save some time.

## Value

A fitted_dlm object.

## Details

This is the main function of the kDGLM package, as it is used to fit all
models.

For the details about the implementation see dos Santos et al. (2024) .

For the details about the methodology see Alves et al. (2024) .

For the details about the Dynamic Linear Models see West and Harrison
(1997) and Petris et al. (2009) .

## See also

auxiliary functions for creating outcomes [`Poisson`](Poisson.md),
[`Multinom`](Multinom.md), [`Normal`](Normal.md), [`Gamma`](Gamma.md)

auxiliary functions for creating structural blocks
[`polynomial_block`](polynomial_block.md),
[`regression_block`](regression_block.md),
[`harmonic_block`](harmonic_block.md), [`TF_block`](tf_block.md)

auxiliary functions for defining priors
[`zero_sum_prior`](zero_sum_prior.md), [`CAR_prior`](CAR_prior.md)

Other auxiliary functions for fitted_dlm objects:
[`coef.fitted_dlm()`](coef.fitted_dlm.md),
[`eval_dlm_norm_const()`](eval_dlm_norm_const.md),
[`forecast.fitted_dlm()`](forecast.fitted_dlm.md),
[`kdglm()`](kdglm.md),
[`simulate.fitted_dlm()`](simulate.fitted_dlm.md),
[`smoothing()`](smoothing.md),
[`update.fitted_dlm()`](update.fitted_dlm.md)

## Examples

``` r
# Poisson case
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = 0.95)
season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)

outcome <- Poisson(lambda = "rate", data = data)

fitted.data <- fit_model(level, season,
  AirPassengers = outcome
)
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     AirPassengers: Poisson
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : -580.2514
#> Interval Score        : 129.48462
#> Mean Abs. Scaled Error:   0.47013
#> ---

plot(fitted.data, plot.pkg = "base")


##################################################################

# Multinomial case
structure <- (
  polynomial_block(p = 1, order = 2, D = 0.95) +
    harmonic_block(p = 1, period = 12, D = 0.975) +
    noise_block(p = 1, R1 = 0.1) +
    regression_block(p = chickenPox$date >= as.Date("2013-09-01"))
  # Vaccine was introduced in September of 2013
) * 4

outcome <- Multinom(p = structure$pred.names, data = chickenPox[, c(2, 3, 4, 6, 5)])
fitted.data <- fit_model(structure, chickenPox = outcome)
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     chickenPox: Multinomial
#> 
#> Static coeficients (smoothed):
#>                  Estimate Std. Error   t value Pr(>|t|)
#> Var.Reg.1         0.39743    0.25059  1.58601    0.113   
#> Var.Reg.2         0.47441    0.26448  1.79376    0.073   
#> Var.Reg.3         0.48811    0.28497  1.71284    0.087   
#> Var.Reg.4        -0.26900    0.23557 -1.14192    0.253   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : -1952.613
#> Interval Score        : 165.55741
#> Mean Abs. Scaled Error:   0.77058
#> ---
plot(fitted.data, plot.pkg = "base")


##################################################################

# Univariate Normal case
structure <- polynomial_block(mu = 1, D = 0.95) +
  polynomial_block(V = 1, D = 0.95)

outcome <- Normal(mu = "mu", V = "V", data = cornWheat$corn.log.return[1:500])
fitted.data <- fit_model(structure, corn = outcome)
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     corn: Normal
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : 1277.54
#> Interval Score        : 0.07564
#> Mean Abs. Scaled Error: 0.76327
#> ---
plot(fitted.data, plot.pkg = "base")


##################################################################

# Bivariate Normal case
structure <- (polynomial_block(mu = 1, D = 0.95) +
  polynomial_block(V = 1, D = 0.95)) * 2 +
  polynomial_block(rho = 1, D = 0.95)

outcome <- Normal(
  mu = c("mu.1", "mu.2"),
  V = matrix(c("V.1", "rho", "rho", "V.2"), 2, 2),
  data = cornWheat[1:500, c(4, 5)]
)
fitted.data <- fit_model(structure, cornWheat = outcome)
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     cornWheat: Normal
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : 2557.802
#> Interval Score        : 0.07505
#> Mean Abs. Scaled Error: 0.74750
#> ---
plot(fitted.data, plot.pkg = "base")


##################################################################

# Gamma case
structure <- polynomial_block(mu = 1, D = 0.95)

Y <- (cornWheat$corn.log.return[1:500] - mean(cornWheat$corn.log.return[1:500]))**2
outcome <- Gamma(phi = 0.5, mu = "mu", data = Y)
fitted.data <- fit_model(structure, corn = outcome)
summary(fitted.data)
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     corn: Gamma
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : 3508.409
#> Interval Score        : 0.00197
#> Mean Abs. Scaled Error: 0.93721
#> ---
plot(fitted.data, plot.pkg = "base")


##################################################################
# \donttest{
# Sensitivity analysis
data <- c(AirPassengers)

level <- polynomial_block(rate = 1, order = 2, D = "D.level")
season <- harmonic_block(rate = "sazo.effect", order = 2, period = 12, D = "D.sazo")

outcome <- Poisson(lambda = "rate", data = data)

fit_model(level, season, outcome,
  sazo.effect = c(0, 1),
  D.level = c(seq.int(0.8, 1, l = 11)),
  D.sazo = c(seq.int(0.95, 1, l = 11)),
  condition = "sazo.effect==1 | D.sazo==1"
)
#> [                                                  ] - 0% - ETA - NA[                                                  ] - 0.76% - ETA - 0.19 minutes[=                                                 ] - 1.52% - ETA - 0.21 minutes[=                                                 ] - 2.27% - ETA - 0.18 minutes[==                                                ] - 3.03% - ETA - 0.18 minutes[==                                                ] - 3.79% - ETA - 0.17 minutes[==                                                ] - 4.55% - ETA - 0.16 minutes[===                                               ] - 5.3% - ETA - 0.15 minutes[===                                               ] - 6.06% - ETA - 0.15 minutes[===                                               ] - 6.82% - ETA - 0.15 minutes[====                                              ] - 7.58% - ETA - 0.15 minutes[====                                              ] - 8.33% - ETA - 0.14 minutes[=====                                             ] - 9.09% - ETA - 0.14 minutes[=====                                             ] - 9.85% - ETA - 0.14 minutes[=====                                             ] - 10.61% - ETA - 0.14 minutes[======                                            ] - 11.36% - ETA - 0.14 minutes[======                                            ] - 12.12% - ETA - 0.14 minutes[======                                            ] - 12.88% - ETA - 0.14 minutes[=======                                           ] - 13.64% - ETA - 0.14 minutes[=======                                           ] - 14.39% - ETA - 0.13 minutes[========                                          ] - 15.15% - ETA - 0.13 minutes[========                                          ] - 15.91% - ETA - 0.13 minutes[========                                          ] - 16.67% - ETA - 0.13 minutes[=========                                         ] - 17.42% - ETA - 0.13 minutes[=========                                         ] - 18.18% - ETA - 0.13 minutes[=========                                         ] - 18.94% - ETA - 0.13 minutes[==========                                        ] - 19.7% - ETA - 0.12 minutes[==========                                        ] - 20.45% - ETA - 0.12 minutes[===========                                       ] - 21.21% - ETA - 0.12 minutes[===========                                       ] - 21.97% - ETA - 0.12 minutes[===========                                       ] - 22.73% - ETA - 0.12 minutes[============                                      ] - 23.48% - ETA - 0.12 minutes[============                                      ] - 24.24% - ETA - 0.11 minutes[============                                      ] - 25% - ETA - 0.11 minutes[=============                                     ] - 25.76% - ETA - 0.11 minutes[=============                                     ] - 26.52% - ETA - 0.11 minutes[==============                                    ] - 27.27% - ETA - 0.11 minutes[==============                                    ] - 28.03% - ETA - 0.11 minutes[==============                                    ] - 28.79% - ETA - 0.11 minutes[===============                                   ] - 29.55% - ETA - 0.11 minutes[===============                                   ] - 30.3% - ETA - 0.11 minutes[================                                  ] - 31.06% - ETA - 0.1 minutes[================                                  ] - 31.82% - ETA - 0.1 minutes[================                                  ] - 32.58% - ETA - 0.1 minutes[=================                                 ] - 33.33% - ETA - 0.1 minutes[=================                                 ] - 34.09% - ETA - 0.1 minutes[=================                                 ] - 34.85% - ETA - 0.1 minutes[==================                                ] - 35.61% - ETA - 0.1 minutes[==================                                ] - 36.36% - ETA - 0.1 minutes[===================                               ] - 37.12% - ETA - 0.1 minutes[===================                               ] - 37.88% - ETA - 0.09 minutes[===================                               ] - 38.64% - ETA - 0.09 minutes[====================                              ] - 39.39% - ETA - 0.09 minutes[====================                              ] - 40.15% - ETA - 0.09 minutes[====================                              ] - 40.91% - ETA - 0.09 minutes[=====================                             ] - 41.67% - ETA - 0.09 minutes[=====================                             ] - 42.42% - ETA - 0.09 minutes[======================                            ] - 43.18% - ETA - 0.09 minutes[======================                            ] - 43.94% - ETA - 0.09 minutes[======================                            ] - 44.7% - ETA - 0.08 minutes[=======================                           ] - 45.45% - ETA - 0.08 minutes[=======================                           ] - 46.21% - ETA - 0.08 minutes[=======================                           ] - 46.97% - ETA - 0.08 minutes[========================                          ] - 47.73% - ETA - 0.08 minutes[========================                          ] - 48.48% - ETA - 0.08 minutes[=========================                         ] - 49.24% - ETA - 0.08 minutes[=========================                         ] - 50% - ETA - 0.08 minutes[=========================                         ] - 50.76% - ETA - 0.08 minutes[==========================                        ] - 51.52% - ETA - 0.07 minutes[==========================                        ] - 52.27% - ETA - 0.07 minutes[===========================                       ] - 53.03% - ETA - 0.07 minutes[===========================                       ] - 53.79% - ETA - 0.07 minutes[===========================                       ] - 54.55% - ETA - 0.07 minutes[============================                      ] - 55.3% - ETA - 0.07 minutes[============================                      ] - 56.06% - ETA - 0.07 minutes[============================                      ] - 56.82% - ETA - 0.07 minutes[=============================                     ] - 57.58% - ETA - 0.06 minutes[=============================                     ] - 58.33% - ETA - 0.06 minutes[==============================                    ] - 59.09% - ETA - 0.06 minutes[==============================                    ] - 59.85% - ETA - 0.06 minutes[==============================                    ] - 60.61% - ETA - 0.06 minutes[===============================                   ] - 61.36% - ETA - 0.06 minutes[===============================                   ] - 62.12% - ETA - 0.06 minutes[===============================                   ] - 62.88% - ETA - 0.06 minutes[================================                  ] - 63.64% - ETA - 0.06 minutes[================================                  ] - 64.39% - ETA - 0.05 minutes[=================================                 ] - 65.15% - ETA - 0.05 minutes[=================================                 ] - 65.91% - ETA - 0.05 minutes[=================================                 ] - 66.67% - ETA - 0.05 minutes[==================================                ] - 67.42% - ETA - 0.05 minutes[==================================                ] - 68.18% - ETA - 0.05 minutes[==================================                ] - 68.94% - ETA - 0.05 minutes[===================================               ] - 69.7% - ETA - 0.05 minutes[===================================               ] - 70.45% - ETA - 0.04 minutes[====================================              ] - 71.21% - ETA - 0.04 minutes[====================================              ] - 71.97% - ETA - 0.04 minutes[====================================              ] - 72.73% - ETA - 0.04 minutes[=====================================             ] - 73.48% - ETA - 0.04 minutes[=====================================             ] - 74.24% - ETA - 0.04 minutes[======================================            ] - 75% - ETA - 0.04 minutes[======================================            ] - 75.76% - ETA - 0.04 minutes[======================================            ] - 76.52% - ETA - 0.04 minutes[=======================================           ] - 77.27% - ETA - 0.03 minutes[=======================================           ] - 78.03% - ETA - 0.03 minutes[=======================================           ] - 78.79% - ETA - 0.03 minutes[========================================          ] - 79.55% - ETA - 0.03 minutes[========================================          ] - 80.3% - ETA - 0.03 minutes[=========================================         ] - 81.06% - ETA - 0.03 minutes[=========================================         ] - 81.82% - ETA - 0.03 minutes[=========================================         ] - 82.58% - ETA - 0.03 minutes[==========================================        ] - 83.33% - ETA - 0.03 minutes[==========================================        ] - 84.09% - ETA - 0.02 minutes[==========================================        ] - 84.85% - ETA - 0.02 minutes[===========================================       ] - 85.61% - ETA - 0.02 minutes[===========================================       ] - 86.36% - ETA - 0.02 minutes[============================================      ] - 87.12% - ETA - 0.02 minutes[============================================      ] - 87.88% - ETA - 0.02 minutes[============================================      ] - 88.64% - ETA - 0.02 minutes[=============================================     ] - 89.39% - ETA - 0.02 minutes[=============================================     ] - 90.15% - ETA - 0.01 minutes[=============================================     ] - 90.91% - ETA - 0.01 minutes[==============================================    ] - 91.67% - ETA - 0.01 minutes[==============================================    ] - 92.42% - ETA - 0.01 minutes[===============================================   ] - 93.18% - ETA - 0.01 minutes[===============================================   ] - 93.94% - ETA - 0.01 minutes[===============================================   ] - 94.7% - ETA - 0.01 minutes[================================================  ] - 95.45% - ETA - 0.01 minutes[================================================  ] - 96.21% - ETA - 0.01 minutes[================================================  ] - 96.97% - ETA - 0 minutes[================================================= ] - 97.73% - ETA - 0 minutes[================================================= ] - 98.48% - ETA - 0 minutes[==================================================] - 99.24% - ETA - 0 minutes
#>     sazo.effect D.level D.sazo  log.like      mase interval.score
#> 240           1    0.98  1.000 -565.4798 0.4109207       104.9385
#> 218           1    0.98  0.995 -565.6162 0.4098146       104.2846
#> 196           1    0.98  0.990 -566.1562 0.4114235       105.2615
#> 174           1    0.98  0.985 -567.0578 0.4145855       106.3615
#> 152           1    0.98  0.980 -568.2681 0.4187701       107.8692
#> Fitted DGLM with 1 outcomes.
#> 
#> distributions:
#>     Series.1: Poisson
#> 
#> ---
#> No static coeficients.
#> ---
#> See the coef.fitted_dlm for the coeficients with temporal dynamic.
#> 
#> One-step-ahead prediction
#> Log-likelihood        : -565.4798
#> Interval Score        : 104.93846
#> Mean Abs. Scaled Error:   0.41092
#> ---
# }
```
