#' Fitting kDGLM models
#'
#' Fit a model given its structure and the observed data. This function can be used for any supported family (see vignette).
#'
#' @param ... dlm_block object: The structural blocks of the model. All block must be completely defined.
#' @param outcomes List: The observed data. It should contain objects of the class dlm_distr.
#' @param pred_cred Numeric: A number between 0 and 1 (not included) indicating the credibility interval for predictions. If not within the valid range of values, predictions are not made.
#' @param smooth Bool: A flag indicating if the smoothed distribution for the latent variables should be calculated.
#' @param p_monit numeric (optional): The prior probability of changes in the latent space variables that are not part of it's dynamic.
#' @param c_monit numeric (optional, if p_monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
#'
#' @return A fitted_dlm object. Contains:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x t.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x t.
#'    \item mts Matrix: The smoothed mean of the latent variables for each time. Dimensions are n x t. Only available if smooth=TRUE.
#'    \item Cts Array: A 3D-array containing the smoothed covariance matrix of the latent variable for each time. Dimensions are n x n x t. Only available if smooth=TRUE.
#'    \item at Matrix: The one-step-ahead mean of the latent variables at each time. Dimensions are n x t.
#'    \item Rt Array: A 3D-array containing the one-step-ahead covariance matrix for latent variables at each time. Dimensions are n x n x t.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are k x t.
#'    \item Qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are k x k x t.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item G_labs Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values) when there is no monitoring. When monitoring for abnormalities, the value in times where abnormalities were detected is increased.
#'    \item h Matrix: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). Its dimension should be n x t, where t is the length of the series and n is the number of latent states.
#'    \item H Array: The same as the argument (same values).
#'    \item W Array: A 3D-array containing the effective covariance matrix of the noise for each time, i.e., considering both H and D. Its dimension should be the same as H and D.
#'    \item outcome List: The same as the argument outcome (same values).
#'    \item pred_names Vector: The names of the linear predictors.
#'    \item var_names List: A list containing names and indexes for latent variables.
#'    \item a1 Matrix: The prior mean for time 1. Dimensions are n x 1.
#'    \item R1 Matrix: The prior covariance matrix for time 1. Dimensions are n x n.
#'    \item smooth Bool: The same as the argument (same value).
#'    \item pred_cred Numeric: The same as the argument (same value).
#'    \item t Numeric: The time range for which the model has been fitted.
#'    \item structure dlm_block: The structure of the model. It's equivalent to block_superpos(...), but also taking into consideration the outcome length.
#' }
#' @export
#'
#' @examples
#' library(kDGLM)
#'
#' # Poisson case
#' T <- 200
#' w <- (T / 40) * 2 * pi
#' data <- rpois(T, exp(sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = 0.95)
#' season <- harmonic_block(rate = 1, period = 40, D = 0.98)
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' ##################################################################
#'
#' # Multinomial case
#' T <- 200
#' y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
#' y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
#' y3 <- rpois(T, exp(5))
#'
#' y <- cbind(y1, y2, y3)
#'
#' level <- polynomial_block(p1 = 1) + polynomial_block(p2 = 1)
#' season <- harmonic_block(p2 = 1, period = 12)
#' outcome <- Multinom(p = c("p1", "p2"), outcome = y)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' ##################################################################
#'
#' # Normal case
#' T <- 200
#' mu <- rnorm(T, 0, 0.1)
#' data <- rnorm(T, cumsum(mu))
#'
#' level <- polynomial_block(
#'   mu = 1,
#'   D = 0.95
#' )
#' variance <- polynomial_block(
#'   sigma2 = 1
#' )
#'
#' # Known variance
#' outcome <- Normal(mu = "mu", Sigma = 1, outcome = data)
#'
#' fitted_data <- fit_model(level, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' # Unknown variance
#' outcome <- Normal(mu = "mu", Sigma = "sigma2", outcome = data)
#'
#' fitted_data <- fit_model(level, variance, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' ##################################################################
#'
#' # Gamma case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' mu <- exp(sin(w * 1:T / T) + 2)
#' data <- matrix(rgamma(T, phi, phi / mu), T, 1)
#'
#' level <- polynomial_block(mu = 1, D = 0.95)
#' season <- harmonic_block(mu = 1, period = 40, D = 0.98)
#' scale <- polynomial_block(phi = 1, D = 1)
#'
#' # Known shape
#' outcome <- Gamma(phi = phi, mu = "mu", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' # DO NOT RUN
#' # # Unknown shape
#' # outcome <- Gamma(phi = "phi", mu = "mu", outcome = data)
#' #
#' # fitted_data <- fit_model(level, season, scale, outcomes = outcome)
#' # summary(fitted_data)
#' #
#' # show_fit(fitted_data, lag=-1)$plot
#'
#' @details
#'
#' This is the main function of the kDGLM package, as it is used to fit all models.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about the methodology see  \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For the details about the Dynamic Linear Models see  \insertCite{WestHarr-DLM;textual}{kDGLM} and \insertCite{Petris-DLM;textual}{kDGLM}.
#'
#' @seealso auxiliary functions for creating outcomes \code{\link{Poisson}}, \code{\link{Multinom}}, \code{\link{Normal}}, \code{\link{Gamma}}, \code{\link{Dirichlet}}
#' @seealso auxiliary functions for creating structural blocks \code{\link{polynomial_block}}, \code{\link{harmonic_block}}, \code{\link{AR_block}}
#' @seealso auxiliary function for choosing hyper parameters \code{\link{search_model}}.
#' @family {auxiliary functions for fitted_dlm objects}
fit_model <- function(..., outcomes, pred_cred = 0.95, smooth = TRUE, p_monit = NA, c_monit = 1) {
  pred_cred <- if (is.numeric(pred_cred)) {
    if (pred_cred > 0 & pred_cred < 1) {
      pred_cred
    } else {
      NA
    }
  } else {
    NA
  }
  if (inherits(outcomes, "dlm_distr")) {
    outcomes <- list(outcomes)
  }

  if (is.null(names(outcomes))) {
    names(outcomes) <- paste0("Serie_", 1:length(outcomes))
  } else if (any(names(outcomes) == "")) {
    val_r <- sum(names(outcomes) == "")
    names(outcomes)[names(outcomes) == ""] <- paste0("Serie_", 1:val_r)
  }

  t <- sapply(outcomes, function(x) {
    x$t
  })
  if (min(t) != max(t)) {
    stop(paste0("Error: outcomes does not have the same time length."))
  }
  t <- max(t)

  structure <- block_superpos(...)
  if (structure$status == "undefined") {
    stop("Error: One or more hiper parameter are undefined. Did you meant to use the search_model function?")
  }

  block_names <- names(structure$var_names)

  for (name in unique(block_names)) {
    count_name <- sum(block_names == name)
    if (count_name > 1) {
      len_char <- floor(log10(count_name)) + 1
      block_names[block_names == name] <- paste0(name, "_", formatC(1:count_name, width = len_char, flag = "0"))
    }
  }
  names(structure$var_names) <- block_names

  coef_names <- rep(NA, structure$n)
  for (name in names(structure$var_names)) {
    coef_names[structure$var_names[[name]]] <- paste0(name, "_", names(structure$var_names[[name]]))
  }
  structure$var_labels <- coef_names

  if (structure$t == 1) {
    structure$t <- t
    structure$G <- array(structure$G, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$G))
    structure$D <- array(structure$D, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$D))
    structure$h <- matrix(structure$h, structure$n, structure$t)
    structure$H <- array(structure$H, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$H))
    structure$FF <- array(structure$FF, c(structure$n, structure$k, structure$t), dimnames = dimnames(structure$FF))
    structure$FF_labs <- matrix(structure$FF_labs, structure$n, structure$k)
  }
  structure$G[, , 1] <- diag(structure$n)
  structure$D[, , 1] <- 1
  structure$H[, , 1] <- 0
  if (t != structure$t) {
    stop(paste0("Error: outcome does not have the same time length as structure: got ", t, " from outcome, expected ", structure$t))
  }

  pred_names_out <- unique(structure$FF_labs)
  pred_names_out <- pred_names_out[pred_names_out != "const"]
  for (outcome in outcomes) {
    pred_names_out <- c(pred_names_out, outcome$pred_names)
  }
  pred_names_out <- unique(pred_names_out)
  if (any(!(pred_names_out %in% structure$pred_names))) {
    stop("Error: One or more linear predictor in outcomes are not present in the model structure.")
  }
  if (any(!(structure$pred_names %in% pred_names_out))) {
    warning("One or more linear predictor in the model structure are not used in the outcomes.")
  }

  model <- analytic_filter(
    outcomes = outcomes,
    a1 = structure$a1,
    R1 = structure$R1,
    FF = structure$FF,
    FF_labs = structure$FF_labs,
    G = structure$G,
    G_labs = structure$G_labs,
    D = structure$D,
    h = structure$h,
    H = structure$H,
    p_monit = p_monit,
    c_monit = c_monit,
    monitoring = structure$monitoring
  )
  if (smooth) {
    model <- smoothing(model)
  }
  for (outcome_name in names(model$outcomes)) {
    outcome <- model$outcomes[[outcome_name]]

    prediction <- outcome$calc_pred(outcome$conj_prior_param, if (is.na(pred_cred)) {
      NULL
    } else {
      outcome$outcome
    }, pred_cred = pred_cred, parms = outcome$parms)

    model$outcomes[[outcome_name]]$pred <- prediction$pred
    model$outcomes[[outcome_name]]$var.pred <- prediction$var.pred
    model$outcomes[[outcome_name]]$icl.pred <- prediction$icl.pred
    model$outcomes[[outcome_name]]$icu.pred <- prediction$icu.pred
    model$outcomes[[outcome_name]]$log.like <- prediction$log.like
  }
  model$a1 <- structure$a1
  model$R1 <- structure$R1
  model$var_names <- structure$var_names
  model$var_labels <- structure$var_labels
  model$pred_cred <- pred_cred
  model$t <- t
  model$k <- structure$k
  model$n <- structure$n
  model$structure <- structure
  class(model) <- "fitted_dlm"

  return(model)
}

#' Auxiliary function for model smoothing
#'
#' @param model A fitted_dlm object.
#'
#' @return A fitted_dlm object with smoothed means (mts) and covariance matrix (Cts) for each observation.
#' @export
smoothing <- function(model) {
  if (!model$smooth) {
    smoothed <- generic_smoother(model$mt, model$Ct, model$at, model$Rt, model$G, model$G_labs)
    model$mts <- smoothed$mts
    model$Cts <- smoothed$Cts
    model$smooth <- TRUE
  } else {
    warning("Model already smoothed.")
  }
  return(model)
}

#' forecast
#'
#' Makes predictions for t times ahead using the fitted model.
#'
#' @param model fitted_dlm: The fitted model to be use for predictions.
#' @param t Numeric: Time window for prediction.
#' @param outcome List (optional): A named list containing the observed values in the prediction window. Note that the names in the list must be the same as the names passed during the fitting process.
#' @param offset Matrix or scalar (optional): offset for predictions. Should have dimensions r x t, where r is the number of outcomes of the model. If offset is not specified, the last value observed by the model will be used.
#' @param FF Array (optional): A 3D-array containing the regression matrix for each time. Its dimension should be n x k x t, where n is the number of latent variables, k is the number of linear predictors in the model. If not specified, the last value given to the model will be used.
#' @param G Array (optional): A 3D-array containing the evolution matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables. If not specified, the last value given to the model will be used.
#' @param D Array (optional): A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and T is the time series length. If not specified, the last value given to the model will be used in the first step, and 1 will be use thereafter.
#' @param h Matrix (optional): A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). Its dimension should be n x t, where t is the length of the series and n is the number of latent states.
#' @param H Array (optional): A 3D-array containing the covariance matrix of the noise for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length. If not specified, 0 will be used.
#' @param plot Bool or String: A flag indicating if a plot should be produced. Should be one of FALSE, TRUE, 'base', 'ggplot2' or 'plotly'.
#' @param pred_cred Numeric: The credibility level for the I.C. intervals.
#'
#' @return A list containing:
#' \itemize{
#'    \item data Data.frame: A data frame contain the mean, variance and credibility intervals for the outcomes, including both the observed data and the predictions for future observations.
#'    \item forecast Data.frame: Same as data, but restricted to predictions for future observations.
#'    \item outcomes List: A named list containing predictions for each outcome. Each element of this list is a list containing predictions (mean, variance and credibility intervals), the distribution of the linear predictor for the parameter of the observational model and the parameters of the predictive distribution (if available).
#'    \item mt Matrix: A matrix with the values for the latent states at each time. Dimensions are n x t, where n is the number of latent states
#'    \item Ct Array: A 3D-array with the covariance of the latent states at each time. Dimensions are n x n x t, where n is the number of linear predictors.
#'    \item ft Matrix: A matrix with the values for the linear predictors at each time. Dimensions are n x t, where n is the number of linear predictors
#'    \item Qt Array: A 3D-array with the covariance of the linear predictors at each time. Dimensions are n x n x t, where n is the number of linear predictors.
#'    \item plot (if so chosen): A plotly or ggplot object.
#' }
#' @importFrom Rfast transpose data.frame.to_matrix
#' @import graphics
#' @import grDevices
#' @export
#'
#' @examples
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = 0.95)
#' season <- harmonic_block(rate = 1, period = 40, D = 0.98)
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#'
#' # forecast(fitted_data, 20)$plot
#' # Or
#' forecast_dlm(fitted_data, 20)$plot
#'
#' @family {auxiliary functions for fitted_dlm objects}
forecast_dlm <- function(model, t = 1, outcome = NULL, offset = NULL, FF = NULL, G = NULL, D = NULL, h = NULL, H = NULL, plot = ifelse(requireNamespace("plotly", quietly = TRUE), "plotly", ifelse(requireNamespace("ggplot2", quietly = TRUE), "ggplot", "base")), pred_cred = 0.95) {
  if (plot == TRUE) {
    plot <- ifelse(requireNamespace("plotly", quietly = TRUE), "plotly", ifelse(requireNamespace("ggplot2", quietly = TRUE), "ggplot", "base"))
  }

  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  k <- dim(model$FF)[2]
  r <- sum(sapply(model$outcomes, function(x) {
    x$r
  }))
  pred.names <- model$pred_names
  pred <- matrix(NA, r, t)
  var.pred <- array(NA, c(r, r, t))
  icl.pred <- matrix(NA, r, t)
  icu.pred <- matrix(NA, r, t)
  log.like <- rep(NA, t)

  #### Consistency check ####
  if (length(dim(FF)) > 3) {
    stop(paste0("Error: FF should have at most 3 dimensions. Got ", length(dim(FF)), "."))
  }
  if (length(dim(G)) > 3) {
    stop(paste0("Error: G should have at most 3 dimensions. Got ", length(dim(G)), "."))
  }
  if (length(dim(D)) > 3) {
    stop(paste0("Error: D should have at most 3 dimensions. Got ", length(dim(D)), "."))
  }
  if (length(dim(H)) > 3) {
    stop(paste0("Error: H should have at most 3 dimensions. Got ", length(dim(H)), "."))
  }
  if (length(dim(offset)) > 2) {
    stop(paste0("Error: D should have at most 2 dimensions. Got ", length(dim(offset)), "."))
  }

  G_labs <- model$G_labs
  if (is.null(G)) {
    G <- array(model$G[, , t_last], c(n, n, t))
  }
  if (is.null(h)) {
    h <- matrix(model$h[, t_last], n, t)
  }
  FF_labs <- model$FF_labs
  if (is.null(FF)) {
    FF <- array(model$FF[, , t_last], c(n, k, t))
  }
  if (is.null(D)) {
    D <- array(1, c(n, n, t))
  } else {
    if (all(dim(D) == 1)) {
      D <- array(D, c(n, n, t))
    } else {
      if (length(dim(D)) == 2 | (length(dim(D)) == 3 & dim(D)[3] == 1)) {
        D <- array(D, c(n, n, t))
      }
    }
  }
  if (is.null(H)) {
    H <- array(model$H[, , t_last], c(n, n, t))
    # H[, , -1] <- 0
  } else {
    if (all(dim(H) == 1)) {
      H <- array(diag(n) * H, c(n, n, t))
    } else {
      if (length(dim(H)) == 2 | (length(dim(H)) == 3 & dim(H)[3] == 1)) {
        H <- array(H, c(n, n, t))
        H[, , -1] <- 0
      }
    }
  }

  if (any(dim(G) != c(n, n, t))) {
    stop(paste0("Error: G has wrong dimesions. Expected ", paste(n, n, t, sep = "x"), ". Got ", paste(dim(G), colapse = "x")))
  }
  if (any(dim(FF) != c(n, k, t))) {
    stop(paste0("Error: FF has wrong dimesions. Expected ", paste(n, k, t, sep = "x"), ". Got ", paste(dim(FF), colapse = "x")))
  }
  if (any(dim(D) != c(n, n, t))) {
    stop(paste0("Error: D has wrong dimesions. Expected ", paste(n, n, t, sep = "x"), ". Got ", paste(dim(D), colapse = "x")))
  }
  if (any(dim(H) != c(n, n, t))) {
    stop(paste0("Error: H has wrong dimesions. Expected ", paste(n, n, t, sep = "x"), ". Got ", paste(dim(H), colapse = "x")))
  }
  #####
  D <- ifelse(D == 0, 1, D)

  a1 <- model$mt[, t_last]
  R1 <- model$Ct[, , t_last]

  m1 <- matrix(NA, n, t)
  C1 <- array(NA, c(n, n, t))
  f1 <- matrix(NA, k, t)
  Q1 <- array(NA, c(k, k, t))

  last_m <- a1
  last_C <- R1

  R0 <- one_step_evolve(
    last_m, last_C,
    G[, , 1] |> matrix(n, n), G_labs,
    D[, , 1]**0, h[, 1], H[, , 1] * 0
  )$Rt

  outcome_forecast <- list()
  for (outcome_name in names(model$outcomes)) {
    r_i <- model$outcomes[[outcome_name]]$r
    outcome_forecast[[outcome_name]]$conj_param <- matrix(NA, t, length(model$outcomes[[outcome_name]]$conj_prior_param)) |> as.data.frame()
    names(outcome_forecast[[outcome_name]]$conj_param) <- names(model$outcomes[[outcome_name]]$conj_prior_param)

    outcome_forecast[[outcome_name]]$ft <- matrix(NA, model$outcomes[[outcome_name]]$k, t)
    outcome_forecast[[outcome_name]]$Qt <- array(NA, c(model$outcomes[[outcome_name]]$k, model$outcomes[[outcome_name]]$k, t))

    if (!is.null(outcome)) {
      outcome_forecast[[outcome_name]]$outcome <- outcome[[outcome_name]] |> matrix(t, r_i)
    } else {
      outcome_forecast[[outcome_name]]$outcome <- model$outcomes[[outcome_name]]$outcome[t_last, ] |> matrix(t, r_i)
    }
    if (!is.null(offset)) {
      outcome_forecast[[outcome_name]]$offset <- offset[[outcome_name]] |> matrix(t, r_i)
    } else {
      outcome_forecast[[outcome_name]]$offset <- model$outcomes[[outcome_name]]$offset[t_last, ] |> matrix(t, r_i)
    }
  }


  for (t_i in c(1:t)) {
    W_i <- as.matrix((1 - D[, , t_i]) * R0 / D[, , t_i])

    next_step <- one_step_evolve(last_m, last_C, G[, , t_i] |> matrix(n, n), G_labs, D[, , t_i]**0, h[, t_i], H[, , t_i] + W_i)
    last_m <- next_step$at
    last_C <- next_step$Rt

    m1[, t_i] <- last_m
    C1[, , t_i] <- last_C

    lin_pred <- calc_lin_pred(last_m, last_C, FF[, , t_i] |> matrix(n, k, dimnames = list(NULL, pred.names)), FF_labs)
    f1[, t_i] <- lin_pred$ft
    Q1[, , t_i] <- lin_pred$Qt
    for (outcome_name in names(model$outcomes)) {
      model_i <- model$outcomes[[outcome_name]]
      pred_index <- match(model_i$pred_names, model$pred_names)
      lin_pred_offset <- model_i$apply_offset(
        lin_pred$ft[pred_index, , drop = FALSE],
        lin_pred$Qt[pred_index, pred_index, drop = FALSE],
        outcome_forecast[[outcome_name]]$offset[t_i, ]
      )
      if (model_i$convert_canom_flag) {
        ft_canom <- model_i$convert_mat_canom %*% lin_pred_offset$ft
        Qt_canom <- model_i$convert_mat_canom %*% lin_pred_offset$Qt %*% transpose(model_i$convert_mat_canom)
      } else {
        ft_canom <- lin_pred_offset$ft
        Qt_canom <- lin_pred_offset$Qt
      }

      conj_prior <- model_i$conj_prior(ft_canom, Qt_canom, parms = model_i$parms)
      outcome_forecast[[outcome_name]]$conj_param[t_i, ] <- conj_prior

      outcome_forecast[[outcome_name]]$ft[, t_i] <- lin_pred_offset$ft
      outcome_forecast[[outcome_name]]$Qt[, , t_i] <- lin_pred_offset$Qt
    }
  }
  if (is.null(outcome)) {
    outcome_forecast[[outcome_name]]$outcome <- NULL
    outcome_forecast[[outcome_name]]$offset <- NULL
  }

  r_acum <- 0
  out_names <- rep(NA, r)
  output <- matrix(NA, t, r)
  for (outcome_name in names(model$outcomes)) {
    model_i <- model$outcomes[[outcome_name]]
    r_cur <- model_i$r

    prediction <- model_i$calc_pred(outcome_forecast[[outcome_name]]$conj_param,
      outcome_forecast[[outcome_name]]$outcome,
      pred_cred = pred_cred,
      parms = model_i$parms
    )


    outcome_forecast[[outcome_name]]$pred <- prediction$pred
    outcome_forecast[[outcome_name]]$var.pred <- prediction$var.pred
    outcome_forecast[[outcome_name]]$icl.pred <- prediction$icl.pred
    outcome_forecast[[outcome_name]]$icu.pred <- prediction$icu.pred

    pred[(r_acum + 1):(r_acum + r_cur), ] <- prediction$pred
    var.pred[(r_acum + 1):(r_acum + r_cur), (r_acum + 1):(r_acum + r_cur), ] <- prediction$var.pred
    icl.pred[(r_acum + 1):(r_acum + r_cur), ] <- prediction$icl.pred
    icu.pred[(r_acum + 1):(r_acum + r_cur), ] <- prediction$icu.pred



    if (r_cur > 1) {
      out_names[(r_acum + 1):(r_acum + r_cur)] <- paste0(outcome_name, "_", 1:r_cur)
    } else {
      out_names[(r_acum + 1):(r_acum + r_cur)] <- outcome_name
    }
    if (!is.null(outcome)) {
      output[, (r_acum + 1):(r_acum + r_cur)] <- outcome_forecast[[outcome_name]]$outcome
    }

    r_acum <- r_acum + r_cur
  }

  data <- data.frame(
    Time = 1:t + t_last,
    Serie = as.factor(c(sapply(out_names, function(x) {
      rep(x, t)
    }))),
    Observation = c(output),
    Prediction = c(t(pred)),
    C.I.lower = c(t(icl.pred)),
    C.I.upper = c(t(icu.pred))
  )

  data_past <- eval_past(model, lag = -1, pred_cred = pred_cred)$data
  plot_data <- rbind(
    cbind(data_past, type = "Fit"),
    cbind(data, type = "Forecast")
  )

  return_list <- list(data = plot_data, forecast = data, outcomes = outcome_forecast, mt = m1, Ct = C1, ft = f1, Qt = Q1)

  if (plot != FALSE) {
    obs_na_rm <- data_past$Observation[!is.na(data_past$Observation)]
    max_value <- evaluate_max(obs_na_rm - min(obs_na_rm))[[3]] + min(obs_na_rm)
    min_value <- -evaluate_max(-(obs_na_rm - max(obs_na_rm)))[[3]] + max(obs_na_rm)
    plot_data$shape_point <- ifelse(plot_data$type == "Forecast", "Future", "Obs.")
    plot_data$group_ribbon <- paste0(plot_data$Serie, plot_data$type)
    series_names <- levels(plot_data$Serie)
    n_series <- length(series_names)
    colors <- rainbow(n_series, s = 0.6)
    points <- paste0(colors, "55")
    fills <- paste0(colors, "33")
    names(colors) <- names(points) <- names(fills) <- series_names
    if (plot == "base" | !requireNamespace("ggplot2", quietly = TRUE)) {
      if (plot != "base") {
        warning("The ggplot2 package is required for ggplot2 and plotly plots and was not found. Falling back to R base plot functions.")
      }
      cur_height <- dev.size("cm")[2]
      count_spaces <- ceiling(n_series / 4)
      font_cm <- 0.35

      config <- par()
      layout(
        mat = matrix(c(1, 2, 3), 3, 1),
        heights = c(
          cur_height - font_cm * count_spaces - 1 - 0.75,
          0.75,
          font_cm * count_spaces + 1
        )
      )

      par(mar = c(4.1, 4.1, 4.1, 2.1), cex = 1)
      plot(0, 0, type = "n", xlim = c(1, t + t_last), ylim = c(min_value, max_value), ylab = "$y_t$", xlab = "Time")
      for (serie in series_names) {
        plot_serie <- plot_data[plot_data$Serie == serie, ]
        points(plot_serie$Time[1:t_last], plot_serie$Observation[1:t_last],
          col = points[[serie]],
          pch = 16
        )
        points(plot_serie$Time[t_last:(t_last + t)], plot_serie$Observation[t_last:(t_last + t)],
          col = points[[serie]],
          pch = 17
        )
        lines(plot_serie$Time[1:t_last], plot_serie$Prediction[1:t_last], col = colors[[serie]])
        lines(plot_serie$Time[t_last:(t_last + t)], plot_serie$Prediction[t_last:(t_last + t)], col = colors[[serie]], lty = 2)
        base_ribbon(plot_serie$Time, plot_serie$C.I.lower, plot_serie$C.I.upper,
          col = fills[[serie]], lty = 0
        )
      }

      par(mar = c(0, 0, 0, 0), cex = 1)
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
      legend(
        legend = c("Fit", "Forecast", "Obs.", "Future"),
        col = c("black", "black", "black", "black"),
        lty = c(1, 2, 0, 0),
        pch = c(0, 0, 16, 17),
        pt.cex = c(0, 0, 1, 1),
        fill = c("gray", "gray", "white", "white"),
        border = "#ffffff00",
        seg.len = 0.6,
        x = 0.5, y = 1, xjust = 0.5, inset = 0, cex = 0.75, bty = "n", horiz = TRUE
      )
      par(mar = c(0, 0, 0, 0), cex = 1)
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
      legend(
        legend = series_names,
        col = colors,
        lty = rep(0, n_series),
        pch = rep(22, n_series),
        pt.cex = rep(2, n_series),
        pt.bg = colors,
        x = 0.5, xjust = 0.5, y = 1, inset = 0, cex = 0.75, bty = "n",
        ncol = min(4, ceiling(n_series / count_spaces))
      )
      par(mar = config$mar)
    } else {
      plt.obj <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "Time")) +
        ggplot2::geom_line(ggplot2::aes_string(y = "Prediction", linetype = "type", color = "Serie")) +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "C.I.lower", ymax = "C.I.upper", fill = "Serie", group = "group_ribbon"), alpha = 0.25, color = NA) +
        ggplot2::geom_point(ggplot2::aes_string(y = "Observation", shape = "shape_point", color = "Serie")) +
        ggplot2::scale_fill_hue("", na.value = NA) +
        ggplot2::scale_color_hue("", na.value = NA) +
        ggplot2::scale_linetype_manual("", values = c("solid", "dashed")) +
        ggplot2::scale_shape_manual("", values = c(17, 16)) +
        ggplot2::scale_x_continuous("Time") +
        ggplot2::scale_y_continuous("$y_t$") +
        ggplot2::theme_bw() +
        ggplot2::coord_cartesian(ylim = c(min_value, max_value))
      if (plot == "plotly") {
        if (!requireNamespace("plotly", quietly = TRUE)) {
          warning("The plotly package is required for plotly plots.")
        } else {
          series_names <- unique(plot_data$Serie)
          n_series <- length(series_names)
          plt.obj <- plotly::ggplotly(plt.obj)
          for (i in (1:n_series) - 1) {
            plt.obj$x$data[[2 * i + 1]]$legendgroup <-
              plt.obj$x$data[[2 * i + 1]]$name <-
              plt.obj$x$data[[2 * i + 2]]$legendgroup <-
              plt.obj$x$data[[2 * i + 2]]$name <-
              plt.obj$x$data[[i + 2 * n_series + 1]]$legendgroup <-
              plt.obj$x$data[[i + 2 * n_series + 1]]$name <-
              paste0(series_names[i + 1], ": fitted values")

            plt.obj$x$data[[2 * i + 1 + 3 * n_series]]$legendgroup <-
              plt.obj$x$data[[2 * i + 1 + 3 * n_series]]$name <-
              plt.obj$x$data[[2 * i + 2 + 3 * n_series]]$legendgroup <-
              plt.obj$x$data[[2 * i + 2 + 3 * n_series]]$name <- paste0(series_names[i + 1], ": observations")

            plt.obj$x$data[[2 * i + 1]]$showlegend <-
              plt.obj$x$data[[i + 1 + 2 * n_series]]$showlegend <-
              plt.obj$x$data[[2 * i + 1 + 3 * n_series]]$showlegend <- FALSE
          }
        }
      }
      return_list$plot <- plt.obj
    }
  }

  return(return_list)
}

#' eval_past
#'
#' Evaluates the predictive values for the observed values used to fit the model and it's latent variables.
#' Predictions can be made with smoothed values or with filtered values with a time offset.
#'
#' @param model fitted_dlm: The fitted model to be use for evaluation.
#' @param T vector: A vector of positive integers indicating the time index from which to extract predictions. Default is to extract all observed times in model.
#' @param lag positive integer: The relative offset for forecast. Values for time t will be calculated based on the filtered values of time t-h. If lag is negative, then the smoothed distribuition will be used.
#' @param pred_cred Numeric: The credibility level for the I.C. intervals.
#' @param eval_pred Bool: A flag indicating if the predictions should be calculated.
#'
#' @return A list containing:
#' \itemize{
#'    \item data data.frame: A table with the model evaluated at each observed time.
#'    \item mt Matrix: The mean of the latent variables at each time. Dimensions are n x t, where t is the size of T.
#'    \item Ct Array: A 3D-array containing the covariance matrix of the latent variable at each time. Dimensions are n x n x t.
#'    \item ft Matrix: The mean of the linear predictor at each time. Dimensions are k x t.
#'    \item Qt Array: A 3D-array containing the covariance matrix for the linear predictor at each time. Dimensions are k x k x t.
#' }
#' @importFrom Rfast transpose data.frame.to_matrix
#' @export
#'
#' @examples
#' T <- 200
#' w <- ((T + 20) / 40) * 2 * pi
#' y1 <- matrix(rpois((T + 20), 20 * (sin(w * 1:(T + 20) / (T + 20)) + 2)), (T + 20), 1)
#' y2 <- matrix(rpois((T + 20), 1:(T + 20) / (T + 20) + 1), (T + 20), 1)
#' y3 <- matrix(rpois((T + 20), 6), (T + 20), 1)
#' y <- cbind(y1, y2, y3)
#' y_pred <- y[T:(T + 20), ]
#'
#' y <- y[1:T, ]
#'
#' level <- polynomial_block(p1 = 1) + polynomial_block(p2 = 1)
#' season_2 <- harmonic_block(p2 = 1, period = 20)
#'
#' outcome <- Multinom(p = c("p1", "p2"), outcome = y)
#'
#' fitted_data <- fit_model(level, season_2, outcomes = outcome)
#'
#' past <- eval_past(fitted_data, lag = -1)
#'
#' @family {auxiliary functions for fitted_dlm objects}
eval_past <- function(model, T = 1:t_last, lag = -1, pred_cred = 0.95, eval_pred = TRUE) {
  if (round(lag) != lag) {
    stop(paste0("Error: lag should be a integer. Got ", lag, "."))
  }
  pred.names <- model$pred_names
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  k <- dim(model$FF)[2]
  r <- sum(sapply(model$outcomes, function(x) {
    x$r
  }))

  if (min(T) < 1) {
    warning("Cannot evaluate parameters before time 1.")
    T <- T[T >= 1]
  }
  if (max(T) > t_last) {
    warning("Cannot evaluate parameters after last observation.")
    T <- T[T <= t_last]
  }
  init_ref <- min(T) - lag
  init_t <- min(T)
  final_t <- max(T)
  len_t <- final_t - init_t + 1

  FF <- model$FF
  FF_labs <- model$FF_labs

  mt.pred <- matrix(NA, n, len_t)
  Ct.pred <- array(NA, c(n, n, len_t))
  ft.pred <- matrix(NA, k, len_t)
  Qt.pred <- array(NA, c(k, k, len_t))

  pred <- matrix(NA, r, len_t)
  var.pred <- array(NA, c(r, r, len_t))
  icl.pred <- matrix(NA, r, len_t)
  icu.pred <- matrix(NA, r, len_t)
  log.like <- rep(0, len_t)
  mae <- rep(0, len_t)
  rae <- rep(0, len_t)
  mse <- rep(0, len_t)
  interval.score <- rep(0, len_t)

  D <- model$D
  D_inv <- 1 / D
  h <- model$h
  H <- model$H
  G <- model$G

  if (lag < 0) {
    if (!model$smooth) {
      model <- smoothing(model)
    }

    lag <- 0
    ref_mt <- model$mts
    ref_Ct <- model$Cts
  } else {
    ref_mt <- model$mt
    ref_Ct <- model$Ct
  }
  G_labs <- model$G_labs

  for (i in c(init_t:final_t)) {
    if (i <= lag) {
      mt <- model$a1
      Ct <- model$R1
      lag_i <- i - 1
    } else if (i <= t_last) {
      mt <- ref_mt[, (i - lag):(i - lag)]
      Ct <- ref_Ct[, , (i - lag):(i - lag)]
      lag_i <- lag
    } else {
      mt <- ref_mt[, t_last]
      Ct <- ref_Ct[, , t_last]
      lag_i <- lag + i - t_last
    }
    next_step <- list("at" = mt, "Rt" = Ct)
    if (lag_i >= 1) {
      for (t in c(1:lag_i)) {
        next_step <- one_step_evolve(next_step$at, next_step$Rt, G[, , i - t + 1], G_labs, D_inv[, , i - t + 1]**(t == 1), h[, i - t + 1], H[, , i - t + 1])
      }
    }
    lin_pred <- calc_lin_pred(next_step$at |> matrix(n, 1), next_step$Rt, FF[, , i] |> matrix(n, k, dimnames = list(NULL, pred.names)), FF_labs)

    mt.pred[, i - init_t + 1] <- next_step$at
    Ct.pred[, , i - init_t + 1] <- next_step$Rt
    ft.pred[, i - init_t + 1] <- lin_pred$ft
    Qt.pred[, , i - init_t + 1] <- lin_pred$Qt

    r_acum <- 0
    if (eval_pred) {
      for (outcome in model$outcomes) {
        r_cur <- outcome$r

        pred_index <- match(outcome$pred_names, model$pred_names)
        r_i <- length(pred_index)

        cur_step <- outcome$apply_offset(lin_pred$ft[pred_index, , drop = FALSE], lin_pred$Qt[pred_index, pred_index, drop = FALSE], outcome$offset[i, ])

        if (outcome$convert_canom_flag) {
          ft_canom <- outcome$convert_mat_canom %*% cur_step$ft
          Qt_canom <- outcome$convert_mat_canom %*% cur_step$Qt %*% transpose(outcome$convert_mat_canom)
        } else {
          ft_canom <- cur_step$ft
          Qt_canom <- cur_step$Qt
        }


        conj_distr <- outcome$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)
        prediction <- outcome$calc_pred(conj_distr,
          if (i > t_last) {
            NULL
          } else {
            outcome$outcome[i, , drop = FALSE]
          },
          pred_cred,
          parms = outcome$parms
        )
        out_ref <- t(outcome$outcome)[, i, drop = FALSE]

        pred[(r_acum + 1):(r_acum + r_cur), i - init_t + 1] <- prediction$pred
        var.pred[(r_acum + 1):(r_acum + r_cur), (r_acum + 1):(r_acum + r_cur), i - init_t + 1] <- prediction$var.pred
        icl.pred[(r_acum + 1):(r_acum + r_cur), i - init_t + 1] <- prediction$icl.pred
        icu.pred[(r_acum + 1):(r_acum + r_cur), i - init_t + 1] <- prediction$icu.pred
        log.like[i - init_t + 1] <- log.like[i - init_t + 1] + sum(prediction$log.like, na.rm = TRUE)
        mae[i - init_t + 1] <- mae[i - init_t + 1] + sum(abs(out_ref - prediction$pred))
        rae[i - init_t + 1] <- rae[i - init_t + 1] + sum(abs(out_ref - prediction$pred) / ifelse(out_ref == 0, 1, out_ref))
        mse[i - init_t + 1] <- mse[i - init_t + 1] + sum((out_ref - prediction$pred)**2)
        interval.score[i - init_t + 1] <- interval.score[i - init_t + 1] +
          sum((prediction$icu.pred - prediction$icl.pred) +
            2 / (1 - pred_cred) * (prediction$icl.pred - out_ref) * (out_ref < prediction$icl.pred) +
            2 / (1 - pred_cred) * (out_ref - prediction$icu.pred) * (out_ref > prediction$icu.pred))
        r_acum <- r_acum + r_cur
      }
    }
  }

  r_acum <- 0
  out_names <- rep(NA, r)
  output <- matrix(NA, len_t, r)
  for (outcome_name in names(model$outcomes)) {
    r_cur <- model$outcomes[[outcome_name]]$r
    char_len <- floor(log10(r_cur)) + 1
    if (r_cur > 1) {
      out_names[(r_acum + 1):(r_acum + r_cur)] <- paste0(outcome_name, "_", formatC(1:r_cur, width = char_len, flag = "0"))
    } else {
      out_names[(r_acum + 1):(r_acum + r_cur)] <- outcome_name
    }
    output[1:min(final_t - init_t + 1, t_last - init_t + 1), (r_acum + 1):(r_acum + r_cur)] <- model$outcomes[[outcome_name]]$outcome[init_t:min(final_t, t_last), ]


    r_acum <- r_acum + r_cur
  }

  data <- data.frame(
    Time = init_t:final_t,
    Serie = as.factor(c(sapply(out_names, function(x) {
      rep(x, final_t - init_t + 1)
    }))),
    Observation = c(output),
    Prediction = c(t(pred)),
    C.I.lower = c(t(icl.pred)),
    C.I.upper = c(t(icu.pred))
  )

  return(list(
    data = data[data$Time %in% T, ],
    mt = mt.pred[, init_t:final_t %in% T, drop = FALSE],
    Ct = Ct.pred[, , init_t:final_t %in% T, drop = FALSE],
    ft = ft.pred[, init_t:final_t %in% T, drop = FALSE],
    Qt = Qt.pred[, , init_t:final_t %in% T, drop = FALSE],
    log.like = log.like[init_t:final_t %in% T, drop = FALSE],
    mae = mae[init_t:final_t %in% T, drop = FALSE],
    rae = rae[init_t:final_t %in% T, drop = FALSE],
    mse = mse[init_t:final_t %in% T, drop = FALSE],
    interval.score = interval.score[init_t:final_t %in% T, drop = FALSE]
  ))
}

#' dlm_sampling
#'
#' This is function draws samples from the latent states using the backward sampling algorithm. See \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 15, for details.
#'
#' @param model fitted_dlm: A fitted model from which to sample.
#' @param sample_size integer: The number of samples to draw.
#' @param filtered_distr bool: A flag indicating if the samples should be dawned from the filtered (if TRUE) or smoothed distribution. The default is FALSE.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Array: An array containing the samples of latent states. Dimensions are n x t x sample_size, where n is the number of latent variables in the model and T is the number of observed values.
#'    \item ft Array: An array containing the samples of linear predictors. Dimensions are k x t x sample_size, where k is the number of linear predictors in the model and t is the number of observed values.
#'    \item param List: A named list containing, for each model outcome, an array with the samples of the parameters of the observational model. Each array will have dimensions l x t x sample_size, where l is the number of parameters in the observational model and t is the number of observed values.
#' }
#'
#' @importFrom Rfast transpose
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#'
#' T <- 200
#' mu <- rnorm(T, 0, 0.1)
#' data <- rnorm(T, cumsum(mu))
#'
#' level <- polynomial_block(
#'   mu = 1,
#'   D = 0.95
#' )
#' variance <- polynomial_block(
#'   sigma2 = 1
#' )
#' # Unknown variance
#' outcome <- Normal(mu = "mu", Sigma = "sigma2", outcome = data)
#'
#' fitted_data <- fit_model(level, variance, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' sample <- dlm_sampling(fitted_data, 2000)
#'
#' @family {auxiliary functions for fitted_dlm objects}
dlm_sampling <- function(model, sample_size, filtered_distr = FALSE) {
  G <- model$G
  G_labs <- model$G_labs
  at <- model$at
  mt <- model$mt
  FF <- model$FF
  FF_labs <- model$FF_labs
  T_len <- dim(mt)[2]
  n <- dim(mt)[1]
  k <- dim(FF)[2]
  sum(sapply(model$outcomes, function(x) {
    x$r
  }))
  outcomes <- list()
  for (outcome_name in names(model$outcomes)) {
    outcomes[[outcome_name]] <- list(
      inv_link = model$outcomes[[outcome_name]]$inv_link_function,
      apply_offset = model$outcomes[[outcome_name]]$apply_offset,
      offset = model$outcomes[[outcome_name]]$offset,
      l = model$outcomes[[outcome_name]]$l,
      r = model$outcomes[[outcome_name]]$r,
      convert_mat_canom = model$outcomes[[outcome_name]]$convert_mat_canom,
      convert_canom_flag = model$outcomes[[outcome_name]]$convert_canom_flag
    )
  }

  mt_sample <- rnorm(n * T_len * sample_size) |> array(c(n, T_len, sample_size))
  ft_sample <- array(NA, c(k, T_len, sample_size))

  Ct_chol <- var_decomp(model$Ct[, , T_len])
  # mt_sample_i <- transpose(Ct_chol) %*% mt_sample[, T_len, ] + mt[, T_len]
  mt_sample_i <- crossprod(Ct_chol, mt_sample[, T_len, ]) + mt[, T_len]
  FF_step <- FF[, , T_len]
  if (any(is.na(FF_step))) {
    ft_sample_i <- sapply(1:sample_size,
      function(j) {
        calc_lin_pred(mt_sample_i[, j], Ct_chol * 0, FF_step, FF_labs)$ft
      },
      simplify = "matrix"
    )
  } else {
    ft_sample_i <- crossprod(FF_step, mt_sample_i)
  }

  pred_names <- colnames(FF)
  for (outcome_name in names(outcomes)) {
    pred_index <- match(model$outcome[[outcome_name]]$pred_names, pred_names)
    r_i <- length(pred_index)

    outcomes[[outcome_name]]$r_i <- r_i
    outcomes[[outcome_name]]$pred_index <- pred_index

    offset_step <- outcomes[[outcome_name]]$offset[T_len, ]
    if (any(is.na(offset_step))) {
      offset_step <- 1
    }

    if (outcomes[[outcome_name]]$convert_canom_flag) {
      ft_sample_i_canom <- outcomes[[outcome_name]]$convert_mat_canom %*% ft_sample_i
    } else {
      ft_sample_i_canom <- ft_sample_i
    }

    ft_sample_i_out <- outcomes[[outcome_name]]$apply_offset(
      ft_sample_i_canom[outcomes[[outcome_name]]$pred_index, , drop = FALSE],
      diag(k) * 0,
      offset_step
    )$ft
    param_sample_i <- outcomes[[outcome_name]]$inv_link(ft_sample_i_out)

    outcomes[[outcome_name]]$param_sample <- array(NA, c(outcomes[[outcome_name]]$l, T_len, sample_size))
    outcomes[[outcome_name]]$param_sample[, T_len, ] <- param_sample_i
  }
  mt_sample[, T_len, ] <- mt_sample_i
  ft_sample[, T_len, ] <- ft_sample_i

  for (t in (T_len - 1):1) {
    Rt <- model$Rt[, , t + 1]
    Ct <- model$Ct[, , t]

    if (filtered_distr) {
      mts <- mt[, t]
      Cts <- Ct
      Ct_chol <- var_decomp(Cts)
      mt_sample_i <- crossprod(Ct_chol, mt_sample[, t, ]) + mts
    } else {
      mt_now <- mt[, t]
      Ct_now <- Ct

      G_ref <- calc_current_G(mt_now, Ct_now, G[, , t + 1], G_labs)$G
      simple_Rt_inv <- Ct %*% crossprod(G_ref, ginv(Rt))
      simple_Rt_inv_t <- transpose(simple_Rt_inv)
      # simple_Rt_inv_t <- solve(Rt, G_ref %*% Ct)
      # simple_Rt_inv <- transpose(simple_Rt_inv_t)


      mts <- mt[, t] + simple_Rt_inv %*% (mt_sample_i - at[, t + 1])
      Cts <- Ct - simple_Rt_inv %*% Rt %*% simple_Rt_inv_t
      # Cts <- Ct - crossprod(simple_Rt_inv_t, Rt) %*% simple_Rt_inv_t
      Ct_chol <- var_decomp(Cts)
      mt_sample_i <- crossprod(Ct_chol, mt_sample[, t, ]) + mts
    }

    FF_step <- FF[, , t]
    if (any(is.na(FF_step))) {
      ft_sample_i <- sapply(1:sample_size,
        function(j) {
          calc_lin_pred(mt_sample_i[, j], Ct_chol * 0, FF_step, FF_labs)$ft
        },
        simplify = "matrix"
      )
    } else {
      ft_sample_i <- crossprod(FF_step, mt_sample_i)
    }

    for (outcome_name in names(outcomes)) {
      offset_step <- outcomes[[outcome_name]]$offset[t, ]
      if (any(is.na(offset_step))) {
        offset_step <- 1
      }
      if (outcomes[[outcome_name]]$convert_canom_flag) {
        ft_sample_i_canom <- outcomes[[outcome_name]]$convert_mat_canom %*% ft_sample_i
      } else {
        ft_sample_i_canom <- ft_sample_i
      }

      ft_sample_i_out <- outcomes[[outcome_name]]$apply_offset(
        ft_sample_i_canom[outcomes[[outcome_name]]$pred_index, , drop = FALSE],
        diag(k) * 0,
        offset_step
      )$ft
      param_sample_i <- outcomes[[outcome_name]]$inv_link(ft_sample_i_out)
      outcomes[[outcome_name]]$param_sample[, t, ] <- param_sample_i
    }
    mt_sample[, t, ] <- mt_sample_i
    ft_sample[, t, ] <- ft_sample_i
  }
  return(list(
    "mt" = mt_sample,
    "ft" = ft_sample,
    "param" = lapply(outcomes, function(x) {
      x$param_sample
    })
  ))
}


#' Search the best hyper parameters fo a kDGLM model
#'
#' @param ... dlm_block object: The structural blocks of the model. At least one block must be "undefined". See polynomial_block for more details.
#' @param outcomes List: The observed data. It should contain objects of the class dlm_distr.
#' @param search_grid List: A named list containing the possible values of each undefined hyper parameter.
#' @param condition Character: A character defining which combinations of undefined hyper parameter should be tested. See example for details.
#' @param metric String: The name of the metric to use for model selection. One of log-likelihood for the one-step-ahead prediction ("log.like"), Mean Absolute Error ("mae"), Relative Absolute Error ("rae"), Mean Squared Error ("mse") or Interval Score ("interval.score") \insertCite{interval_score}{kDGLM}.
#' @param smooth Bool: A flag indicating if the smoothed distribution for the latent variables should be calculated.
#' @param lag Integer: The number of steps ahead used for the prediction when calculating the metrics. If lag<0, predictions are made using the smoothed distribuition of the latent variables.
#' @param pred_cred Numeric: A number between 0 and 1 (not included) indicating the credibility interval for predictions. If not within the valid range of values, 0.95 will be used.
#' @param metric_cutoff Integer: The number of observations to ignore when calculating the metrics. Default is 1/10 of the number of observations (rounded down).
#' @param p_monit numeric (optional): The prior probability of changes in the latent space variables that are not part of it's dynamic.
#' @param c_monit numeric (optional, if p_monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
#'
#' @return A searched_dlm object containing the following values:
#' \itemize{
#'    \item search_data Data.frame: A table contain the metrics of each combination of hyper parameter tested and ordered from best combination to worst.
#'    \item structure dlm_block: The initial structure passed by the user, by with the undefined hyper parameters set to the best values found.
#'    \item model fitted_dlm: A model fitted with the best combination of hyper parameters.
#' }
#'
#' @export
#'
#' @examples
#'
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = "D_level")
#' season <- harmonic_block(rate = "sazo_effect", period = 40, D = "D_sazo")
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' search_model(level, season,
#'   outcomes = outcome,
#'   search_grid = list(
#'     sazo_effect = c(0, 1),
#'     D_level = c(seq(0.8, 1, 0.05)),
#'     D_sazo = c(seq(0.95, 1, 0.01))
#'   ),
#'   condition = "sazo_effect==1 | D_sazo==1"
#' )
#'
#' @details
#'
#' This is an auxiliary function to search for the best combination of hyper parameters among the possible variations specified by the user.
#' This function simple evaluates all possible combinations and chooses the one that results in the best fitting.
#'
#' @seealso \code{\link{fit_model}}
#' @references
#'    \insertAllCited{}
search_model <- function(..., outcomes, search_grid, condition = "TRUE", metric = "log.like", smooth = TRUE, lag = 1, pred_cred = 0.95, metric_cutoff = NA, p_monit = NA, c_monit = 1) {
  if (is.null(pred_cred) | is.na(pred_cred)) {
    pred_cred <- 0.95
  } else {
    if (pred_cred <= 0 | pred_cred >= 1) {
      warning("Invalid value for credibility. Using 95% of credibility.")
      pred_cred <- 0.95
    }
  }
  structure <- block_superpos(...)
  if (structure$status == "defined") {
    stop("Error: There is no hiper parameter to select. Did you forgot to label the hiper parameters?")
  }

  if (any(names(search_grid) %in% c("const", "constrained", "free", "kl"))) {
    stop("Error: Invalid label for hyper parameter. Cannot use 'const', 'constrained', 'free' or 'kl'. Chose another label.")
  }

  ref_strucuture <- structure

  if (any(names(search_grid) %in% ref_strucuture$pred_names)) {
    stop(paste0(
      "Error: Ambiguous label for hyper parameter. Cannot have a hyper parameter with the same name of a linear predictor\n
         The user passed the following hyper parameters: ", paste0(names(search_grid), collapse = ", "), ".\n",
      "The model has the the following linear predictors: ", paste0(ref_strucuture$pred_names, collapse = ", "), "."
    ))
  }

  var_length <- length(search_grid)
  search_data <- do.call(expand.grid, search_grid)
  search_data <- search_data[eval(parse(text = condition), envir = search_data), ]

  search_data$log.like <- NA
  search_data$mae <- NA
  search_data$rae <- NA
  search_data$mse <- NA
  search_data$interval.score <- NA
  vals_names <- names(search_grid)

  dim_FF <- dim(structure$FF)
  dimnames_FF <- dimnames(structure$FF)

  dim_FF_labs <- dim(structure$FF_labs)

  dim_D <- dim(structure$D)
  dim_H <- dim(structure$H)
  dim_G <- dim(structure$G)

  dim_R1 <- dim(structure$R1)

  vals_nums <- dim(search_data)[1]
  vals_size <- dim(search_data)[2]
  best_structure <- NULL
  best_model <- NULL
  best_metric <- Inf
  init <- Sys.time()
  for (i in 1:vals_nums) {
    time_past <- Sys.time() - init
    raw_perc <- i / vals_nums
    perc <- round(100 * raw_perc, 2)
    n_bar <- round(50 * raw_perc)
    cat(paste0(
      "\r[", paste0(rep("=", n_bar), collapse = ""),
      paste0(rep(" ", 50 - n_bar), collapse = ""), "] - ",
      perc,
      "% - ETA - ",
      round(as.numeric((1 - raw_perc) * time_past / raw_perc, units = "mins"), 2),
      " minutes                 "
    ))
    cur_param <- search_data[i, ]
    structure <- ref_strucuture
    for (name in vals_names) {
      structure$FF[array(structure$FF_labs == name, dim_FF)] <- cur_param[[name]]
      structure$FF_labs[structure$FF_labs == name] <- "const"
      structure$D[structure$D == name] <- cur_param[[name]]
      structure$H[structure$H == name] <- cur_param[[name]]
      structure$a1[structure$a1 == name] <- cur_param[[name]]
      structure$R1[structure$R1 == name] <- cur_param[[name]]
      structure$G[structure$G == name] <- cur_param[[name]]
    }

    structure$FF <- array(as.numeric(structure$FF), dim_FF, dimnames = dimnames_FF)
    structure$D <- array(as.numeric(structure$D), dim_D)
    structure$H <- array(as.numeric(structure$H), dim_H)
    structure$a1 <- as.numeric(structure$a1)
    structure$R1 <- matrix(as.numeric(structure$R1), dim_R1[1], dim_R1[2])
    structure$G <- array(as.numeric(structure$G), dim_G)
    structure$status <- check.block.status(structure)

    if (structure$status == "undefined") {
      stop("Error: not all unkown hiper parameter have values. Check the search grid to make sure every unkown hiper parameter has a range of values.")
    }

    if (any(if.na(structure$D, 0) < 0 | if.na(structure$D, 0) > 1)) {
      stop(paste0("Error: invalid value for D. Expected a real number between 0 and 1, got: ", paste(structure$D[if.na(structure$D, 0) < 0 | if.na(structure$D, 0) > 1], collapse = ", "), "."))
    }
    if (any(if.na(structure$H, 0) < 0)) {
      stop(paste0("Error: invalid value for H. Expected a non negative number, got: ", paste(structure$H[if.na(structure$H, 0) < 0], collapse = ", "), "."))
    }
    if (any(if.na(structure$R1, 0) < 0)) {
      stop(paste0("Error: invalid value for R1. Expected a non negative number, got: ", paste(structure$R1[if.na(structure$R1, 0) < 0], collapse = ", "), "."))
    }

    fitted_model <- fit_model(structure, outcomes = outcomes, pred_cred = NA, smooth = FALSE, p_monit = p_monit, c_monit = c_monit)

    T <- fitted_model$t
    if (is.na(metric_cutoff)) {
      metric_cutoff <- floor(T / 10)
    }

    r <- length(fitted_model$outcomes)
    metric <- tolower(metric)
    predictions <- eval_past(fitted_model, T = metric_cutoff:T, lag = lag, pred_cred = pred_cred, eval_pred = TRUE)

    search_data$log.like[i] <- sum(predictions$log.like)
    search_data$mae[i] <- mean(predictions$mae)
    search_data$rae[i] <- mean(predictions$rae)
    search_data$mse[i] <- mean(predictions$mse)
    search_data$interval.score[i] <- sum(predictions$interval.score)

    cur_metric <- ifelse(metric == "log.like", -search_data$log.like[i], search_data[[metric]][i])
    if (cur_metric < best_metric) {
      best_structure <- structure
      best_model <- fitted_model
      best_metric <- cur_metric
    }
  }
  cat("\n")
  if (metric == "log.like") {
    search_data <- search_data[order(-search_data[[metric]]), ]
  } else {
    search_data <- search_data[order(search_data[[metric]]), ]
  }
  if (smooth) {
    best_model <- smoothing(best_model)
  }

  out.vals <- list(
    search.data = search_data,
    structure = best_structure,
    model = best_model
  )
  class(out.vals) <- "searched_dlm"

  return(out.vals)
}
