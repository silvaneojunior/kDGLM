# Fitting one kDGLM models

Fits one model given its structure and the observed data. This function
can be used for any supported family (see vignette).

## Usage

``` r
fit_model_single(
  structure,
  outcomes,
  smooth = TRUE,
  p.monit = NA,
  safe.mode = TRUE
)
```

## Arguments

- structure:

  dlm_block: The structural blocks of the model. All block must be
  completely defined.

- outcomes:

  dlm_distr or list of dlm_distr objects: The model outcomes.

- smooth:

  boolean: A flag indicating if the smoothed distribution for the latent
  states should be calculated.

- p.monit:

  numeric (optional): The prior probability of changes in the latent
  space variables that are not part of its dynamic.

- safe.mode:

  boolean: A flag indicating if consistency check should be performed at
  each time step. Recommended to be left on, but if you know what you
  are doing (i.e., you tested the model and it is safe) and need to fit
  it several times, you can disable the checks to save some time.

## Value

A fitted_dlm object.
