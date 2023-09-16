#' Fitting kDGLM models
#'
#' Fit a model given its structure and the observed data. This function can be used for any supported family (see vignette).
#'
#' @param ... dlm_block or dlm_distr objects: The structural blocks of the model (dlm_block objects), alongside the model outcomes (dlm_distr object). All block must be completely defined.
#' @param pred.cred Numeric: A number between 0 and 1 (not included) indicating the credibility interval for predictions. If not within the valid range of values, predictions are not made.
#' @param smooth Bool: A flag indicating if the smoothed distribution for the latent variables should be calculated.
#' @param p.monit numeric (optional): The prior probability of changes in the latent space variables that are not part of its dynamic.
#' @param c.monit numeric (optional, if p.monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
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
#'    \item G.labs Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values) when there is no monitoring. When monitoring for abnormalities, the value in times where abnormalities were detected is increased.
#'    \item h Matrix: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). Its dimension should be n x t, where t is the length of the series and n is the number of latent states.
#'    \item H Array: The same as the argument (same values).
#'    \item W Array: A 3D-array containing the effective covariance matrix of the noise for each time, i.e., considering both H and D. Its dimension should be the same as H and D.
#'    \item outcomes List: A list containing the outcomes of the model passed in the ... argument.
#'    \item pred.names Vector: The names of the linear predictors.
#'    \item var.names List: A list containing names and indexes for latent variables.
#'    \item a1 Matrix: The prior mean for time 1. Dimensions are n x 1.
#'    \item R1 Matrix: The prior covariance matrix for time 1. Dimensions are n x n.
#'    \item smooth Bool: The same as the argument (same value).
#'    \item pred.cred Numeric: The same as the argument (same value).
#'    \item t Numeric: The time range for which the model has been fitted.
#'    \item n Numeric: The number of latent variables in the model.
#'    \item k Numeric: The number of linear predictors in the model.
#'    \item l Numeric: The total number of parameters for the observational model.
#'    \item r Numeric: The total number number of outcomes of the model.
#'    \item structure dlm_block: The structure of the model. It's equivalent to block_superpos(...), but also taking into consideration the outcome length.
#' }
#' @export
#'
#' @examples
#'
#' # Poisson case
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 0.95)
#' season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)
#'
#' outcome <- Poisson(lambda = "rate", data = data)
#'
#' fitted.data <- fit_model(level, season,
#'   AirPassengers = outcome
#' )
#' summary(fitted.data)
#'
#' plot(fitted.data, plot.pkg = "base")
#'
#' ##################################################################
#'
#' # Multinomial case
#' structure <- (
#'   polynomial_block(p = 1, order = 2, D = 0.95) +
#'     harmonic_block(p = 1, period = 12, D = 0.975) +
#'     noise_block(p=1,R1=0.1)+
#'     regression_block(p = as.Date('2013-09-1')) # Vaccine was introduced in September of 2013
#' ) * 4
#'
#' outcome <- Multinom(p = structure$pred.names, data = chickenPox[, c(2, 3, 4, 6, 5)])
#' fitted.data <- fit_model(structure, outcome)
#' summary(fitted.data)
#' plot(fitted.data, plot.pkg = "base")
#'
#' ##################################################################
#'
#' # Univariate Normal case
#' structure <- polynomial_block(mu = 1, D = 0.95) +
#'   polynomial_block(V = 1, D = 0.95)
#'
#' outcome <- Normal(mu = "mu", V = "V", data = cornWheat$corn.log.return[1:500])
#' fitted.data <- fit_model(structure, outcome)
#' summary(fitted.data)
#' plot(fitted.data, plot.pkg = "base")
#'
#' ##################################################################
#'
#' # Bivariate Normal case
#' structure <- (polynomial_block(mu = 1, D = 0.95) +
#'   polynomial_block(V = 1, D = 0.95)) * 2 +
#'   polynomial_block(rho = 1, D = 0.95)
#'
#' outcome <- Normal(
#'   mu = c("mu.1", "mu.2"),
#'   V = matrix(c("V.1", "rho", "rho", "V.2"), 2, 2),
#'   data = cornWheat[1:500, c(4, 5)]
#' )
#' fitted.data <- fit_model(structure, outcome)
#' summary(fitted.data)
#' plot(fitted.data, plot.pkg = "base")
#'
#' ##################################################################
#'
#' # Gamma case
#' structure <- polynomial_block(mu = 1, D = 0.95)
#'
#' outcome <- Gamma(phi = 0.5, mu = "mu", data = cornWheat$corn.log.return[1:500]**2)
#' fitted.data <- fit_model(structure, outcome)
#' summary(fitted.data)
#' plot(fitted.data, plot.pkg = "base")
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
#' @seealso auxiliary functions for creating structural blocks \code{\link{polynomial_block}}, \code{\link{regression_block}}, \code{\link{harmonic_block}}, \code{\link{AR_block}}
#' @seealso auxiliary function for choosing hyper parameters \code{\link{search_model}}.
#' @family {auxiliary functions for fitted_dlm objects}.
fit_model <- function(..., smooth = TRUE, p.monit = NA, c.monit = 1) {

  extra.args <- list(...)
  structure <- list()
  outcomes <- list()
  out.names <- c()
  for (i in seq_along(extra.args)) {
    arg <- extra.args[[i]]
    arg.name <- if.null(names(extra.args)[i], "")
    if (inherits(arg, "dlm_distr")) {
      out.names <- c(out.names, arg.name)
      outcomes[[length(outcomes) + 1]] <- arg
    } else if (inherits(arg, "dlm_block")) {
      structure[[length(structure) + 1]] <- arg
    } else {
      stop(paste0("Error: Invalid type for ... argument. Expected a dlm_block or dlm_distr object, got ", class(arg), ". Be sure that all arguments are properly named."))
    }
  }
  if (length(outcomes) == 0) {
    stop("Error: No dlm_distr object was passed. Make sure all outcomes are created using the proper functions (see documentation).")
  }
  if (length(structure) == 0) {
    stop("Error: No dlm_block object was passed. Make sure all blocks are created using the proper functions (see documentation).")
  }
  if (any(out.names == "")) {
    out.names[out.names == ""] <- paste0("Series.", 1:length(outcomes))[out.names == ""]
  }
  names(outcomes) <- out.names

  structure <- do.call(block_superpos, structure)
  if (structure$status == "undefined") {
    stop("Error: One or more hiper parameter are undefined. Did you meant to use the search_model function?")
  }

  t <- sapply(outcomes, function(x) {
    x$t
  })
  if (min(t) != max(t)) {
    stop(paste0("Error: outcomes does not have the same time length."))
  }
  t <- max(t)

  block.names <- names(structure$var.names)

  for (name in unique(block.names)) {
    count.name <- sum(block.names == name)
    if (count.name > 1) {
      len.char <- floor(log10(count.name)) + 1
      block.names[block.names == name] <- paste0(name, ".", formatC(1:count.name, width = len.char, flag = "0"))
    }
  }
  names(structure$var.names) <- block.names

  coef.names <- rep(NA, structure$n)
  for (name in names(structure$var.names)) {
    coef.names[structure$var.names[[name]]] <- paste0(name, ".", names(structure$var.names[[name]]))
  }
  structure$var.labels <- coef.names

  if (structure$t == 1) {
    structure$t <- t
    structure$G <- array(structure$G, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$G))
    structure$D <- array(structure$D, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$D))
    structure$h <- matrix(structure$h, structure$n, structure$t)
    structure$H <- array(structure$H, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$H))
    structure$FF <- array(structure$FF, c(structure$n, structure$k, structure$t), dimnames = dimnames(structure$FF))
    structure$FF.labs <- matrix(structure$FF.labs, structure$n, structure$k)
  }
  structure$G[, , 1] <- diag(structure$n)
  structure$D[, , 1] <- 1
  structure$h[, 1] <- 0
  structure$H[, , 1] <- 0
  if (t != structure$t) {
    stop(paste0("Error: outcome does not have the same time length as structure: got ", t, " from outcome, expected ", structure$t))
  }

  pred.names.out <- unique(structure$FF.labs)
  pred.names.out <- pred.names.out[pred.names.out != "const"]
  for (outcome in outcomes) {
    pred.names.out <- c(pred.names.out, outcome$pred.names)
  }
  pred.names.out <- unique(pred.names.out)
  if (any(!(pred.names.out %in% structure$pred.names))) {
    stop("Error: One or more linear predictor in outcomes are not present in the model structure.")
  }
  if (any(!(structure$pred.names %in% pred.names.out))) {
    warning("One or more linear predictor in the model structure are not used in the outcomes.")
  }

  for (intervention in structure$interventions) {
    var.index.inter <- intervention$var.index
    if (!is.null(intervention$FF)) {
      structure$FF[
        var.index.inter,
        structure$pred.names %in% intervention$pred.names,
        intervention$times
      ] <- intervention$FF
    }

    if (!is.null(intervention$D)) {
      structure$D[
        var.index.inter, var.index.inter,
        intervention$times
      ] <- intervention$D
    }

    if (!is.null(intervention$h)) {
      structure$h[
        var.index.inter,
        intervention$times
      ] <- intervention$h
    }

    if (!is.null(intervention$H)) {
      structure$H[
        var.index.inter, var.index.inter,
        intervention$times
      ] <- intervention$H
    }

    if (!is.null(intervention$G)) {
      structure$G[
        var.index.inter, var.index.inter,
        intervention$times
      ] <- intervention$G
    }
  }

  model <- analytic_filter(
    outcomes = outcomes,
    a1 = structure$a1,
    R1 = structure$R1,
    FF = structure$FF,
    FF.labs = structure$FF.labs,
    G = structure$G,
    G.labs = structure$G.labs,
    D = structure$D,
    h = structure$h,
    H = structure$H,
    p.monit = p.monit,
    c.monit = c.monit,
    monitoring = structure$monitoring
  )
  if (smooth) {
    model <- smoothing(model)
  }

  r <- sum(sapply(model$outcomes,function(x){x$r}))
  model$a1 <- structure$a1
  model$R1 <- structure$R1
  model$var.names <- structure$var.names
  model$var.labels <- structure$var.labels
  model$t <- t
  model$n <- structure$n
  model$k <- structure$k
  model$l <- structure$l
  model$r <- r
  model$structure <- structure
  model$period <- structure$period
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
    smoothed <- generic_smoother(model$mt, model$Ct, model$at, model$Rt, model$G, model$G.labs)
    model$mts <- smoothed$mts
    model$Cts <- smoothed$Cts

    # model$=eval_dlm_prior(fitted.data$mts,fitted.data)+
    #   eval_dlm_log_like(fitted.data$mts,fitted.data)+
    #   -eval_dlm_post(fitted.data$mts,fitted.data)

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
#' @param D Array (optional): A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length. If not specified, the last value given to the model will be used in the first step, and 1 will be use thereafter.
#' @param h Matrix (optional): A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). Its dimension should be n x t, where t is the length of the series and n is the number of latent states.
#' @param H Array (optional): A 3D-array containing the covariance matrix of the noise for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length. If not specified, 0 will be used.
#' @param plot Bool or String: A flag indicating if a plot should be produced. Should be one of FALSE, TRUE, 'base', 'ggplot2' or 'plotly'.
#' @param pred.cred Numeric: The credibility level for the I.C. intervals.
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
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 0.95)
#' season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)
#'
#' outcome <- Poisson(lambda = "rate", data)
#'
#' fitted.data <- fit_model(level, season,
#'   AirPassengers = outcome
#' )
#'
#' # forecast(fitted.data, 20)$plot
#' # Or
#' forecast_dlm(fitted.data, 20)$plot
#'
#' @family {auxiliary functions for fitted_dlm objects}
forecast_dlm <- function(model, t = 1,
                         outcome = NULL, offset = NULL,
                         FF = NULL, G = NULL,
                         D = NULL, h = NULL, H = NULL,
                         plot = ifelse(requireNamespace("plotly", quietly = TRUE), "plotly", ifelse(requireNamespace("ggplot2", quietly = TRUE), "ggplot2", "base")),
                         pred.cred = 0.95) {
  if (plot == TRUE) {
    plot <- ifelse(requireNamespace("plotly", quietly = TRUE), "plotly", ifelse(requireNamespace("ggplot2", quietly = TRUE), "ggplot2", "base"))
  }

  n <- model$n
  t_last <- model$t
  k <- model$k
  r <- model$r
  pred.names <- model$pred.names
  pred <- matrix(NA, r, t)
  var.pred <- array(NA, c(r, r, t))
  icl.pred <- matrix(NA, r, t)
  icu.pred <- matrix(NA, r, t)

  time.index=seq_len(t_last)
  time.index.foward=(t_last+1):(t_last + t)

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

  G.labs <- model$G.labs
  if (is.null(G)) {
    G <- array(model$G[, , t_last], c(n, n, t))
  }
  if (is.null(h)) {
    h <- matrix(model$h[, t_last], n, t)
  }
  FF.labs <- model$FF.labs
  if (is.null(FF)) {
    FF <- array(model$FF[, , t_last], c(n, k, t))
  }
  dim.D=dim(D)
  if (is.null(D)) {
    D <- array(1, c(n, n, t))
  } else {
    if (all(dim.D == 1)) {
      D <- array(D, c(n, n, t))
    } else {
      if (length(dim.D) == 2 || (length(dim.D) == 3 && dim.D[3] == 1)) {
        D <- array(D, c(n, n, t))
      }
    }
  }
  dim.H=dim(H)
  if (is.null(H)) {
    H <- array(model$H[, , t_last], c(n, n, t))
    # H[, , -1] <- 0
  } else {
    if (all(dim.H == 1)) {
      H <- array(diag(n) * H, c(n, n, t))
    } else {
      if (length(dim.H) == 2 || (length(dim.H) == 3 && dim.H[3] == 1)) {
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

  last.m <- a1
  last.C <- R1

  R0 <- one_step_evolve(
    last.m, last.C,
    G[, , 1] |> matrix(n, n), G.labs,
    D[, , 1]**0, h[, 1], H[, , 1] * 0
  )$Rt

  outcome.forecast <- list()
  for (outcome.name in names(model$outcomes)) {
    r_i <- model$outcomes[[outcome.name]]$r
    outcome.forecast[[outcome.name]]$conj.param <- matrix(NA, t, length(model$outcomes[[outcome.name]]$param.names)) |> as.data.frame()
    names(outcome.forecast[[outcome.name]]$conj.param) <- model$outcomes[[outcome.name]]$param.names
    row.names(outcome.forecast[[outcome.name]]$conj.param )= time.index.foward

    outcome.forecast[[outcome.name]]$ft <- matrix(NA, model$outcomes[[outcome.name]]$k, t)
    outcome.forecast[[outcome.name]]$Qt <- array(NA, c(model$outcomes[[outcome.name]]$k, model$outcomes[[outcome.name]]$k, t))

    outcome.forecast[[outcome.name]]$data <-
      if (!is.null(outcome)) {
        outcome[[outcome.name]] |> matrix(t, r_i)
      } else {
        model$outcomes[[outcome.name]]$data[t_last, ] |> matrix(t, r_i,byrow=TRUE)
        }
    outcome.forecast[[outcome.name]]$offset <-
      if (!is.null(offset)) {offset[[outcome.name]] |> matrix(t, r_i)
      } else {model$outcomes[[outcome.name]]$offset[t_last, ] |> matrix(t, r_i,byrow=TRUE)}
  }


  for (t_i in seq_len(t)) {
    W_i <- matrix((1 - D[, , t_i]) * R0 / D[, , t_i],n,n)

    next.step <- one_step_evolve(last.m, last.C, G[, , t_i] |> matrix(n, n), G.labs, D[, , t_i]**0, h[, t_i], H[, , t_i] + W_i)
    last.m <- next.step$at
    last.C <- next.step$Rt

    m1[, t_i] <- last.m
    C1[, , t_i] <- last.C

    lin.pred <- calc_lin_pred(last.m, last.C,
                              FF[, , t_i] |> matrix(n, k, dimnames = list(NULL, pred.names)),
                              FF.labs, pred.names)
    f1[, t_i] <- lin.pred$ft
    Q1[, , t_i] <- lin.pred$Qt
    for (outcome.name in names(model$outcomes)) {
      model_i <- model$outcomes[[outcome.name]]
      pred.index <- match(model_i$pred.names, model$pred.names)
      lin.pred.offset <- model_i$apply_offset(
        lin.pred$ft[pred.index, , drop = FALSE],
        lin.pred$Qt[pred.index, pred.index, drop = FALSE],
        outcome.forecast[[outcome.name]]$offset[t_i, ]
      )
      if (model_i$convert.canom.flag) {
        ft.canom <- model_i$convert.mat.canom %*% lin.pred.offset$ft
        Qt.canom <- model_i$convert.mat.canom %*% lin.pred.offset$Qt %*% transpose(model_i$convert.mat.canom)
      } else {
        ft.canom <- lin.pred.offset$ft
        Qt.canom <- lin.pred.offset$Qt
      }

      conj.param <- model_i$conj_distr(ft.canom, Qt.canom, parms = model_i$parms)
      outcome.forecast[[outcome.name]]$conj.param[t_i, ] <- conj.param

      outcome.forecast[[outcome.name]]$ft[, t_i] <- lin.pred.offset$ft
      outcome.forecast[[outcome.name]]$Qt[, , t_i] <- lin.pred.offset$Qt
    }
  }
  # if (is.null(outcome)) {
  #   outcome.forecast[[outcome.name]]$data <- matrix(model$outcomes[[outcome.name]]$data[1,],t,,byrow=TRUE)
  #   outcome.forecast[[outcome.name]]$offset <- NULL
  # }

  r.acum <- 0
  out.names <- rep(NA, r)
  output <- matrix(NA, t, r)
  for (outcome.name in names(model$outcomes)) {
    model_i <- model$outcomes[[outcome.name]]
    r.cur <- model_i$r
    prediction <- model_i$calc_pred(outcome.forecast[[outcome.name]]$conj.param,
      outcome.forecast[[outcome.name]]$data,
      pred.cred = pred.cred,
      parms = model_i$parms
    )


    outcome.forecast[[outcome.name]]$pred <- prediction$pred
    outcome.forecast[[outcome.name]]$var.pred <- prediction$var.pred
    outcome.forecast[[outcome.name]]$icl.pred <- prediction$icl.pred
    outcome.forecast[[outcome.name]]$icu.pred <- prediction$icu.pred

    r.seq=(r.acum + 1):(r.acum + r.cur)

    pred[r.seq, ] <- prediction$pred
    var.pred[r.seq, r.seq, ] <- prediction$var.pred
    icl.pred[r.seq, ] <- prediction$icl.pred
    icu.pred[r.seq, ] <- prediction$icu.pred



    if (r.cur > 1) {
      out.names[r.seq] <- paste0(outcome.name, ".", seq_len(r.cur))
    } else {
      out.names[r.seq] <- outcome.name
    }
    if (!is.null(outcome)) {
      output[, r.seq] <- outcome.forecast[[outcome.name]]$data
    }

    r.acum <- r.acum + r.cur
  }

  data <- data.frame(
    Time = time.index.foward,
    Serie = as.factor(c(sapply(out.names, function(x) {
      rep(x, t)
    }))),
    Observation = c(output),
    Prediction = c(t(pred)),
    C.I.lower = c(t(icl.pred)),
    C.I.upper = c(t(icu.pred))
  )

  data.past <- eval_past(model, lag = -1, pred.cred = pred.cred)$data
  plot.data <- rbind(
    cbind(data.past, type = "Fit"),
    cbind(data, type = "Forecast")
  )

  return.list <- list(data = plot.data, forecast = data, outcomes = outcome.forecast, mt = m1, Ct = C1, ft = f1, Qt = Q1)

  if (plot != FALSE) {
    obs.na.rm <- data.past$Observation[!is.na(data.past$Observation)]
    max.value <- evaluate_max(obs.na.rm - min(obs.na.rm))[[3]] + min(obs.na.rm)
    min.value <- -evaluate_max(-(obs.na.rm - max(obs.na.rm)))[[3]] + max(obs.na.rm)
    plot.data$shape.point <- ifelse(plot.data$type == "Forecast", "Future", "Obs.")
    plot.data$group.ribbon <- paste0(plot.data$Serie, plot.data$type)
    series.names <- levels(plot.data$Serie)
    n.series <- length(series.names)
    colors <- rainbow(n.series, s = 0.8)
    fills <- rainbow(n.series, s = 0.4)
    points <- paste0(colors, "55")
    fills <- paste0(fills, "33")
    names(colors) <- names(points) <- names(fills) <- series.names
    if (plot == "base" || !requireNamespace("ggplot2", quietly = TRUE)) {
      if (plot != "base") {
        warning("The ggplot2 package is required for ggplot2 and plotly plots and was not found. Falling back to R base plot functions.")
      }
      cur.height <- dev.size("cm")[2]
      count.spaces <- ceiling(n.series / 4)
      font.cm <- 0.35

      config <- par()
      layout(
        mat = matrix(c(1, 2, 3), 3, 1),
        heights = c(
          cur.height - font.cm * count.spaces - 1 - 0.75,
          0.75,
          font.cm * count.spaces + 1
        )
      )

      par(mar = c(4.1, 4.1, 4.1, 2.1), cex = 1)
      plot(0, 0, type = "n", xlim = c(1, t + t_last), ylim = c(min.value, max.value), ylab = expression(Y[t]), xlab = "Time")


      for (serie in series.names) {
        plot.serie <- plot.data[plot.data$Serie == serie, ]
        points(plot.serie$Time, plot.serie$Observation,
          col = points[[serie]],
          pch = 16+(seq_len(t+t_last)>t_last)
        )

        lines(plot.serie$Time[time.index], plot.serie$Prediction[time.index], col = colors[[serie]],lty=1)
        lines(plot.serie$Time[time.index.foward], plot.serie$Prediction[time.index.foward], col = colors[[serie]],lty=2)
        base_ribbon(plot.serie$Time, plot.serie$C.I.lower, plot.serie$C.I.upper,
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
        legend = series.names,
        col = colors,
        lty = rep(0, n.series),
        pch = rep(22, n.series),
        pt.cex = rep(2, n.series),
        pt.bg = colors,
        x = 0.5, xjust = 0.5, y = 1, inset = 0, cex = 0.75, bty = "n",
        ncol = min(4, ceiling(n.series / count.spaces))
      )
      par(mar = config$mar)
    } else {
      plt.obj <- ggplot2::ggplot(plot.data, ggplot2::aes_string(x = "Time")) +
        ggplot2::geom_line(ggplot2::aes_string(y = "Prediction", linetype = "type", color = "Serie")) +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "C.I.lower", ymax = "C.I.upper", fill = "Serie", group = "group.ribbon"), alpha = 0.25, color = NA) +
        ggplot2::geom_point(ggplot2::aes_string(y = "Observation", shape = "shape.point", color = "Serie")) +
        ggplot2::scale_fill_hue("", na.value = NA) +
        ggplot2::scale_color_hue("", na.value = NA) +
        ggplot2::scale_linetype_manual("", values = c("solid", "dashed")) +
        ggplot2::scale_shape_manual("", values = c(17, 16)) +
        ggplot2::scale_x_continuous("Time") +
        ggplot2::scale_y_continuous(expression(Y[t])) +
        ggplot2::theme_bw() +
        ggplot2::coord_cartesian(ylim = c(min.value, max.value))
      if (plot == "plotly") {
        if (!requireNamespace("plotly", quietly = TRUE)) {
          warning("The plotly package is required for plotly plots.")
        } else {
          series.names <- unique(plot.data$Serie)
          n.series <- length(series.names)
          plt.obj <- plotly::ggplotly(plt.obj+ggplot2::scale_y_continuous(plotly::TeX('Y_t'))) |>plotly::config(mathjax = 'cdn')
          for (i in (1:n.series) - 1) {
            plt.obj$x$data[[2 * i + 1]]$legendgroup <-
              plt.obj$x$data[[2 * i + 1]]$name <-
              plt.obj$x$data[[2 * i + 2]]$legendgroup <-
              plt.obj$x$data[[2 * i + 2]]$name <-
              plt.obj$x$data[[i + 2 * n.series + 1]]$legendgroup <-
              plt.obj$x$data[[i + 2 * n.series + 1]]$name <-
              paste0(series.names[i + 1], ": fitted values")

            plt.obj$x$data[[2 * i + 1 + 3 * n.series]]$legendgroup <-
              plt.obj$x$data[[2 * i + 1 + 3 * n.series]]$name <-
              plt.obj$x$data[[2 * i + 2 + 3 * n.series]]$legendgroup <-
              plt.obj$x$data[[2 * i + 2 + 3 * n.series]]$name <- paste0(series.names[i + 1], ": observations")

            plt.obj$x$data[[2 * i + 1]]$showlegend <-
              plt.obj$x$data[[i + 1 + 2 * n.series]]$showlegend <-
              plt.obj$x$data[[2 * i + 1 + 3 * n.series]]$showlegend <- FALSE
          }
        }
      }
      return.list$plot <- plt.obj
    }
  }

  return(return.list)
}

#' eval_past
#'
#' Evaluates the predictive values for the observed values used to fit the model and its latent variables.
#' Predictions can be made with smoothed values or with filtered values with a time offset.
#'
#' @param model fitted_dlm: The fitted model to be use for evaluation.
#' @param eval_t vector: A vector of positive integers indicating the time index from which to extract predictions. Default is to extract all observed times in model.
#' @param lag positive integer: The relative offset for forecast. Values for time t will be calculated based on the filtered values of time t-h. If lag is negative, then the smoothed distribution will be used.
#' @param pred.cred Numeric: The credibility level for the I.C. intervals.
#' @param eval.pred Bool: A flag indicating if the predictions should be calculated.
#'
#' @return A list containing:
#' \itemize{
#'    \item data data.frame: A table with the model evaluated at each observed time.
#'    \item mt Matrix: The mean of the latent variables at each time. Dimensions are n x t, where t is the size of eval_t.
#'    \item Ct Array: A 3D-array containing the covariance matrix of the latent variable at each time. Dimensions are n x n x t.
#'    \item ft Matrix: The mean of the linear predictor at each time. Dimensions are k x t.
#'    \item Qt Array: A 3D-array containing the covariance matrix for the linear predictor at each time. Dimensions are k x k x t.
#'    \item log.like, mae, mase, rae, mse, interval.score: The metric value at each time.
#'    \item conj.param list: A list containing, for each outcome, a data.frame with the parameter of the conjugated distribution at each time.
#' }
#' @importFrom Rfast transpose data.frame.to_matrix
#' @export
#'
#' @examples
#' # Poisson case
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 0.95)
#' season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)
#'
#' outcome <- Poisson(lambda = "rate", data = data)
#'
#' fitted.data <- fit_model(level, season,
#'   AirPassengers = outcome
#' )
#'
#' # var.vals <- coef(fitted.data)
#' # Or
#' # var.vals <- coefficients(fitted.data)
#' # Or
#' var.vals <- eval_past(fitted.data)
#'
#' @family {auxiliary functions for fitted_dlm objects}
eval_past <- function(model, eval_t = seq_len(model$t), lag = -1, pred.cred = 0.95, eval.pred = TRUE) {
  if (round(lag) != lag) {
    stop(paste0("Error: lag should be a integer. Got ", lag, "."))
  }
  pred.names <- model$pred.names
  n <- model$n
  t_last <- model$t
  k <- model$k
  r <- model$r

  smoothed.log.like=FALSE
  if (lag < 0) {
    if (!model$smooth) {
      model <- smoothing(model)
    }
    smoothed.log.like=TRUE
    lag <- 0
    ref.mt <- model$mts
    ref.Ct <- model$Cts
  } else {
    ref.mt <- model$mt
    ref.Ct <- model$Ct
  }
  if(lag<=0){
    null.rows.flags=(apply(model$G!=0  | is.na(model$G),3,function(x){rowSums(x,na.rm=TRUE)})==0)
    if(any(null.rows.flags)){
      time_index=colSums(null.rows.flags)>0
      var.W=apply(model$W,3,diag)
      # var.W=matrix(var.W[,model$t],n,model$t)
      ref.mt[null.rows.flags]=0
      for(i in (1:model$t)[time_index]){
        null.rows.flag=null.rows.flags[,i]
        diag(ref.Ct[,,i])[null.rows.flag]=var.W[null.rows.flag,i]
      }
    }
  }

  if (min(eval_t) < 1) {
    warning("Cannot evaluate parameters before time 1.")
    eval_t <- eval_t[eval_t - lag >= 1]
  }
  if (max(eval_t) > t_last) {
    warning("Cannot evaluate parameters after last observation.")
    eval_t <- eval_t[eval_t <= t_last]
  }
  init.t <- min(eval_t)
  final.t <- max(eval_t)
  len.t <- final.t - init.t + 1

  FF <- model$FF
  FF.labs <- model$FF.labs

  mt.pred <- matrix(NA, n, len.t)
  Ct.pred <- array(NA, c(n, n, len.t))
  ft.pred <- matrix(NA, k, len.t)
  Qt.pred <- array(NA, c(k, k, len.t))

  pred <- matrix(NA, r, len.t)
  var.pred <- array(NA, c(r, r, len.t))
  icl.pred <- matrix(NA, r, len.t)
  icu.pred <- matrix(NA, r, len.t)
  log.like <- rep(0, len.t)
  mae <- rep(0, len.t)
  rae <- rep(0, len.t)
  mse <- rep(0, len.t)
  mase <- rep(0, len.t)
  interval.score <- rep(0, len.t)



  conj.param.list=list()
  for (outcome.name in names(model$outcomes)) {
    conj.param.list[[outcome.name]]=matrix(NA,len.t,length(model$outcomes[[outcome.name]]$param.names)) |> as.data.frame()
    names(conj.param.list[[outcome.name]])=model$outcomes[[outcome.name]]$param.names
    row.names(conj.param.list[[outcome.name]])=init.t:final.t
  }

  D <- model$D
  D.inv <- 1 / D
  D.holder <- model$D[,,1]
  h <- model$h
  W <- model$W
  G <- model$G
  G.labs <- model$G.labs

  for (i in c(init.t:final.t)) {
    if (i <= lag) {
      mt <- model$a1
      Ct <- model$R1
      lag_i <- i - 1
    } else if (i <= t_last) {
      mt <- ref.mt[, (i - lag):(i - lag)]
      Ct <- ref.Ct[, , (i - lag):(i - lag)]
      lag_i <- lag
    } else {
      mt <- ref.mt[, t_last]
      Ct <- ref.Ct[, , t_last]
      lag_i <- lag + i - t_last
    }
    next.step <- list("at" = mt, "Rt" = Ct)
    if (lag_i >= 1) {
      for (t in c(lag_i:1)) {
        next.step <- one_step_evolve(next.step$at, next.step$Rt, G[, , i - t + 1], G.labs, D.holder, h[, i - t + 1], W[, , i - t + 1])
      }
    }
    lin.pred <- calc_lin_pred(
      next.step$at |> matrix(n, 1),
      next.step$Rt, FF[, , i] |> matrix(n, k, dimnames = list(NULL, pred.names)),
      FF.labs, pred.names
    )

    mt.pred[, i - init.t + 1] <- next.step$at
    Ct.pred[, , i - init.t + 1] <- next.step$Rt
    ft.pred[, i - init.t + 1] <- lin.pred$ft
    Qt.pred[, , i - init.t + 1] <- lin.pred$Qt

    r.acum <- 0
    if (eval.pred) {
      for (outcome in model$outcomes) {
        r.cur <- outcome$r

        r.seq=(r.acum + 1):(r.acum + r.cur)
        t.index=i - init.t + 1

        pred.index <- match(outcome$pred.names, model$pred.names)
        r_i <- length(pred.index)

        cur.step <- outcome$apply_offset(lin.pred$ft[pred.index, , drop = FALSE], lin.pred$Qt[pred.index, pred.index, drop = FALSE], outcome$offset[i, ])

        if (outcome$convert.canom.flag) {
          ft.canom <- outcome$convert.mat.canom %*% cur.step$ft
          Qt.canom <- outcome$convert.mat.canom %*% cur.step$Qt %*% transpose(outcome$convert.mat.canom)
        } else {
          ft.canom <- cur.step$ft
          Qt.canom <- cur.step$Qt
        }


        conj.param <- outcome$conj_distr(ft.canom, Qt.canom, parms = outcome$parms)
        conj.param.list[[outcome.name]][t.index,]=conj.param
        prediction <- outcome$calc_pred(conj.param,
          if (i > t_last) {
            NULL
          } else {
            outcome$data[i, , drop = FALSE]
          },
          pred.cred,
          parms = outcome$parms
        )
        out.ref <- t(outcome$data)[, i, drop = FALSE]
        if (model$period < model$t) {
          mase.coef <- colMeans(abs(diff(outcome$data, model$period)))
        } else if (model$t > 1) {
          mase.coef <- colMeans(abs(diff(outcome$data, 1)))
        } else {
          mase.coef <- outcome$data[1, ]
        }


        pred[r.seq, t.index] <- prediction$pred
        var.pred[r.seq, r.seq, t.index] <- prediction$var.pred
        icl.pred[r.seq, t.index] <- prediction$icl.pred
        icu.pred[r.seq, t.index] <- prediction$icu.pred
        log.like[t.index] <- log.like[t.index] + sum(prediction$log.like, na.rm = TRUE)
        mae[t.index] <- mae[t.index] + sum(abs(out.ref - prediction$pred))
        rae[t.index] <- rae[t.index] + sum(abs(out.ref - prediction$pred) / ifelse(out.ref == 0, 1, out.ref))
        mse[t.index] <- mse[t.index] + sum((out.ref - prediction$pred)**2)
        interval.score[t.index] <- interval.score[t.index] +
          sum((prediction$icu.pred - prediction$icl.pred) +
            2 / (1 - pred.cred) * (prediction$icl.pred - out.ref) * (out.ref < prediction$icl.pred) +
            2 / (1 - pred.cred) * (out.ref - prediction$icu.pred) * (out.ref > prediction$icu.pred))
        r.acum <- r.acum + r.cur

        mase[t.index] <- mase[t.index] + sum(abs(out.ref - prediction$pred) / mase.coef)
      }
    }
  }

  r.acum <- 0
  out.names <- rep(NA, r)
  output <- matrix(NA, len.t, r)
  for (outcome.name in names(model$outcomes)) {
    r.cur <- model$outcomes[[outcome.name]]$r
    r.seq=(r.acum + 1):(r.acum + r.cur)
    char.len <- floor(log10(r.cur)) + 1
    if (r.cur > 1) {
      out.names[r.seq] <- paste0(outcome.name, ".", formatC(seq_len(r.cur), width = char.len, flag = "0"))
    } else {
      out.names[r.seq] <- outcome.name
    }

    1:min(final.t - init.t + 1, t_last - init.t + 1)
    output[seq_len(min(final.t - init.t + 1,t_last - init.t + 1)), r.seq] <- model$outcomes[[outcome.name]]$data[init.t:min(final.t, t_last), ]


    r.acum <- r.acum + r.cur
  }

  time.index.final=init.t:final.t
  time.flags=time.index.final %in% eval_t

  data <- data.frame(
    Time = time.index.final,
    Serie = as.factor(c(sapply(out.names, function(x) {
      rep(x, final.t - init.t + 1)
    }))),
    Observation = c(output),
    Prediction = c(t(pred)),
    C.I.lower = c(t(icl.pred)),
    C.I.upper = c(t(icu.pred))
  )

  return(list(
    data = data[data$Time %in% eval_t, ],
    mt = mt.pred[, time.flags, drop = FALSE],
    Ct = Ct.pred[, , time.flags, drop = FALSE],
    ft = ft.pred[, time.flags, drop = FALSE],
    Qt = Qt.pred[, , time.flags, drop = FALSE],
    log.like = if(smoothed.log.like){
      # model$mts[,]=ref.mt
      # model$Cts[,,]=ref.Ct
      eval_dlm_prior(model$mts,model)+
        eval_dlm_log_like(model$mts,model)+
        -eval_dlm_post(model$mts,model)
      }else{
      log.like[time.flags, drop = FALSE]
      },
    mae = mae[time.flags, drop = FALSE],
    mase = mase[time.flags, drop = FALSE],
    rae = rae[time.flags, drop = FALSE],
    mse = mse[time.flags, drop = FALSE],
    interval.score = interval.score[time.flags, drop = FALSE],
    conj.param=conj.param.list
  ))
}

#' dlm_sampling
#'
#' This is function draws samples from the latent states using the backward sampling algorithm. See \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 15, for details.
#'
#' @param model fitted_dlm: A fitted model from which to sample.
#' @param sample.size integer: The number of samples to draw.
#' @param filt.distr bool: A flag indicating if the samples should be dawned from the filtered (if TRUE) or smoothed distribution. The default is FALSE.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Array: An array containing the samples of latent states. Dimensions are n x t x sample.size, where n is the number of latent variables in the model and t is the number of observed values.
#'    \item ft Array: An array containing the samples of linear predictors. Dimensions are k x t x sample.size, where k is the number of linear predictors in the model and t is the number of observed values.
#'    \item param List: A named list containing, for each model outcome, an array with the samples of the parameters of the observational model. Each array will have dimensions l x t x sample.size, where l is the number of parameters in the observational model and t is the number of observed values.
#' }
#'
#' @importFrom Rfast transpose
#' @export
#'
#' @examples
#'
#' structure <- polynomial_block(mu = 1, D = 0.95) +
#'   polynomial_block(V = 1, D = 0.95)
#'
#' outcome <- Normal(mu = "mu", V = "V", data = cornWheat$corn.log.return[1:500])
#' fitted.data <- fit_model(structure, outcome)
#'
#' sample <- dlm_sampling(fitted.data, 5000)
#'
#' @family {auxiliary functions for fitted_dlm objects}
dlm_sampling <- function(model, sample.size, filt.distr = FALSE) {
  G <- model$G
  G.labs <- model$G.labs
  at <- model$at
  Rt <- model$Rt
  mt <- model$mt
  Ct <- model$Ct
  FF <- model$FF
  FF.labs <- model$FF.labs
  pred.names=model$pred.names
  t.len <- model$t
  n <- model$n
  k <- model$k
  r <- model$r
  outcomes <- list()
  for (outcome.name in names(model$outcomes)) {
    outcomes[[outcome.name]] <- list(
      inv_link = model$outcomes[[outcome.name]]$inv_link_function,
      apply_offset = model$outcomes[[outcome.name]]$apply_offset,
      offset = model$outcomes[[outcome.name]]$offset,
      l = model$outcomes[[outcome.name]]$l,
      r = model$outcomes[[outcome.name]]$r,
      convert.mat.canom = model$outcomes[[outcome.name]]$convert.mat.canom,
      convert.canom.flag = model$outcomes[[outcome.name]]$convert.canom.flag
    )
  }

  mt.sample <- array(NA, c(n, t.len, sample.size))
  ft.sample <- array(NA, c(k, t.len, sample.size))

  mt.sample_i <- rmvnorm(sample.size, mt[, t.len], Ct[, , t.len])
  FF.step <- FF[, , t.len]
  if (any(is.na(FF.step))) {
    ft.sample_i <- sapply(seq_len(sample.size),
      function(j) {
        calc_lin_pred(mt.sample_i[, j], Ct.chol * 0, FF.step, FF.labs, pred.names)$ft
      },
      simplify = "matrix"
    )
  } else {
    ft.sample_i <- crossprod(FF.step, mt.sample_i)
  }

  pred.names <- colnames(FF)
  for (outcome.name in names(outcomes)) {
    pred.index <- match(model$outcome[[outcome.name]]$pred.names, pred.names)
    r_i <- length(pred.index)

    outcomes[[outcome.name]]$r_i <- r_i
    outcomes[[outcome.name]]$pred.index <- pred.index

    offset.step <- outcomes[[outcome.name]]$offset[t.len, ]
    if (any(is.na(offset.step))) {
      offset.step <- 1
    }

    if (outcomes[[outcome.name]]$convert.canom.flag) {
      ft.sample_i_canom <- outcomes[[outcome.name]]$convert.mat.canom %*% ft.sample_i
    } else {
      ft.sample_i_canom <- ft.sample_i
    }

    ft.sample_i_out <- outcomes[[outcome.name]]$apply_offset(
      ft.sample_i_canom[outcomes[[outcome.name]]$pred.index, , drop = FALSE],
      diag(k) * 0,
      offset.step
    )$ft
    param.sample_i <- outcomes[[outcome.name]]$inv_link(ft.sample_i_out)

    outcomes[[outcome.name]]$param.sample <- array(NA, c(outcomes[[outcome.name]]$l, t.len, sample.size))
    outcomes[[outcome.name]]$param.sample[, t.len, ] <- param.sample_i
  }
  mt.sample[, t.len, ] <- mt.sample_i
  ft.sample[, t.len, ] <- ft.sample_i

  for (t in (t.len - 1):1) {
    Rt.step <- Rt[, , t + 1]
    Ct.step <- Ct[, , t]

    if (filt.distr) {
      mts <- mt[, t]
      Cts <- Ct.step
    } else {
      mt.now <- mt[, t]
      Ct.now <- Ct.step

      G.ref <- calc_current_G(mt.now, Ct.now, G[, , t + 1], G.labs)$G
      simple.Rt.inv <- Ct.step %*% crossprod(G.ref, ginv(Rt.step))
      simple.Rt.inv.t <- transpose(simple.Rt.inv)

      mts <- mt.now + simple.Rt.inv %*% (mt.sample_i - at[, t + 1])
      Cts <- Ct.step - simple.Rt.inv %*% Rt.step %*% simple.Rt.inv.t
    }
    mt.sample_i <- rmvnorm(sample.size, rep(0, n), Cts) + mts

    FF.step <- FF[, , t]
    if (any(is.na(FF.step))) {
      ft.sample_i <- sapply(seq_len(sample.size),
        function(j) {
          calc_lin_pred(mt.sample_i[, j], Ct.chol * 0, FF.step, FF.labs, pred.names)$ft
        },
        simplify = "matrix"
      )
    } else {
      ft.sample_i <- crossprod(FF.step, mt.sample_i)
    }

    for (outcome.name in names(outcomes)) {
      offset.step <- outcomes[[outcome.name]]$offset[t, ]
      if (any(is.na(offset.step))) {
        offset.step <- 1
      }
      if (outcomes[[outcome.name]]$convert.canom.flag) {
        ft.sample_i_canom <- outcomes[[outcome.name]]$convert.mat.canom %*% ft.sample_i
      } else {
        ft.sample_i_canom <- ft.sample_i
      }

      ft.sample_i_out <- outcomes[[outcome.name]]$apply_offset(
        ft.sample_i_canom[outcomes[[outcome.name]]$pred.index, , drop = FALSE],
        diag(k) * 0,
        offset.step
      )$ft
      param.sample_i <- outcomes[[outcome.name]]$inv_link(ft.sample_i_out)
      outcomes[[outcome.name]]$param.sample[, t, ] <- param.sample_i
    }
    mt.sample[, t, ] <- mt.sample_i
    ft.sample[, t, ] <- ft.sample_i
  }
  return(list(
    "mt" = mt.sample,
    "ft" = ft.sample,
    "param" = lapply(outcomes, function(x) {
      x$param.sample
    })
  ))
}


#' Search the best hyper parameters fo a kDGLM model
#'
#' @param ... dlm_block or dlm_distr objects: The structural blocks of the model (dlm_block objects), alongside the model outcomes (dlm_distr object). All block must be completely defined.
#' @param search.grid List: A named list containing the possible values of each undefined hyper parameter.
#' @param condition Character: A character defining which combinations of undefined hyper parameter should be tested. See example for details.
#' @param metric String: The name of the metric to use for model selection. One of log-likelihood for the one-step-ahead prediction ("log.like"), Mean Absolute Error ("mae"), Mean Absolute Scaled Error ("mase") \insertCite{mase}{kDGLM}, Relative Absolute Error ("rae"), Mean Squared Error ("mse") or Interval Score ("interval.score") \insertCite{interval_score}{kDGLM}.
#' @param smooth Bool: A flag indicating if the smoothed distribution for the latent variables should be calculated.
#' @param lag Integer: The number of steps ahead used for the prediction when calculating the metrics. If lag<0, predictions are made using the smoothed distribution of the latent variables.
#' @param pred.cred Numeric: A number between 0 and 1 (not included) indicating the credibility interval for predictions. If not within the valid range of values, 0.95 will be used.
#' @param metric.cutoff Integer: The number of observations to ignore when calculating the metrics. Default is 1/10 of the number of observations (rounded down).
#' @param p.monit numeric (optional): The prior probability of changes in the latent space variables that are not part of its dynamic.
#' @param c.monit numeric (optional, if p.monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
#'
#' @return A searched_dlm object containing the following values:
#' \itemize{
#'    \item search.data Data.frame: A table contain the metrics of each combination of hyper parameter tested and ordered from best combination to worst.
#'    \item structure dlm_block: The initial structure passed by the user, by with the undefined hyper parameters set to the best values found.
#'    \item model fitted_dlm: A model fitted with the best combination of hyper parameters.
#' }
#'
#' @export
#'
#' @examples
#'
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 'D.level')
#' season <- harmonic_block(rate = 'sazo.effect', order = 2, period = 12, D = "D.sazo")
#'
#' outcome <- Poisson(lambda = "rate", data = data)
#'
#' search_model(level, season,
#'   outcomes = outcome,
#'   search.grid = list(
#'     sazo.effect = c(0, 1),
#'     D.level = c(seq.int(0.8, 1, l = 11)),
#'     D.sazo = c(seq.int(0.95, 1, l = 11))
#'   ),
#'   condition = "sazo.effect==1 | D.sazo==1"
#' )
#' @details
#'
#' This is an auxiliary function to search for the best combination of hyper parameters among the possible variations specified by the user.
#' This function simple evaluates all possible combinations and chooses the one that results in the best fitting.
#'
#' @seealso \code{\link{fit_model}}
#' @references
#'    \insertAllCited{}
search_model <- function(..., search.grid, condition = "TRUE", metric = "log.like", smooth = TRUE, lag = 1, pred.cred = 0.95, metric.cutoff = NA, p.monit = NA, c.monit = 1) {
  if (is.null(pred.cred) || is.na(pred.cred)) {
    pred.cred <- 0.95
  } else {
    if (pred.cred <= 0 || pred.cred >= 1) {
      warning("Invalid value for credibility. Using 95% of credibility.")
      pred.cred <- 0.95
    }
  }


  extra.args <- list(...)
  structure <- list()
  outcomes <- list()
  out.names <- c()
  for (i in seq_along(extra.args)) {
    arg <- extra.args[[i]]
    arg.name <- if.null(names(extra.args)[i], "")
    if (inherits(arg, "dlm_distr")) {
      out.names <- c(out.names, arg.name)
      outcomes[[length(outcomes) + 1]] <- arg
    } else if (inherits(arg, "dlm_block")) {
      structure[[length(structure) + 1]] <- arg
    } else {
      stop(paste0("Error: Invalid type for ... argument. Expected a dlm_block or dlm_distr object, got ", class(arg), ". Be sure that all arguments are properly named."))
    }
  }
  if (length(outcomes) == 0) {
    stop("Error: No dlm_distr object was passed. Make sure all outcomes are created using the proper functions (see documentation).")
  }
  if (length(structure) == 0) {
    stop("Error: No dlm_block object was passed. Make sure all blocks are created using the proper functions (see documentation).")
  }
  if (any(out.names == "")) {
    out.names[out.names == ""] <- paste0("Series.", 1:length(outcomes))[out.names == ""]
  }
  names(outcomes) <- out.names

  structure <- do.call(block_superpos, structure)
  if (structure$status == "defined") {
    stop("Error: There is no hiper parameter to select. Did you forgot to label the hiper parameters?")
  }

  if (any(names(search.grid) %in% c("const", "constrained", "free", "kl"))) {
    stop("Error: Invalid label for hyper parameter. Cannot use 'const', 'constrained', 'free' or 'kl'. Chose another label.")
  }

  ref.strucuture <- structure

  if (any(names(search.grid) %in% ref.strucuture$pred.names)) {
    stop(paste0(
      "Error: Ambiguous label for hyper parameter. Cannot have a hyper parameter with the same name of a linear predictor\n
         The user passed the following hyper parameters: ", paste0(names(search.grid), collapse = ", "), ".\n",
      "The model has the the following linear predictors: ", paste0(ref.strucuture$pred.names, collapse = ", "), "."
    ))
  }

  var.length <- length(search.grid)
  search.data <- do.call(expand.grid, search.grid)
  if (condition != "TRUE") {
    search.data <- search.data[eval(parse(text = condition), envir = search.data), ]
  }
  if (dim(search.data)[1] == 0) {
    stop("Error: No model to test. Verify if the search.grid is properly specified and if a valid value for condition argument is being passed.")
  }

  search.data$log.like <- NA
  search.data$mae <- NA
  search.data$rae <- NA
  search.data$mse <- NA
  search.data$mase <- NA
  search.data$interval.score <- NA
  vals.names <- names(search.grid)

  dimnames.FF <- dimnames(structure$FF)

  vals.nums <- dim(search.data)[1]
  vals.size <- dim(search.data)[2]
  best.structure <- NULL
  best.model <- NULL
  best.metric <- Inf
  init <- Sys.time()
  for (i in seq_len(vals.nums)) {
    time.past <- Sys.time() - init
    raw.perc <- i / vals.nums
    perc <- round(100 * raw.perc, 2)
    n.bar <- round(50 * raw.perc)
    cat(paste0(
      "\r[", paste0(rep("=", n.bar), collapse = ""),
      paste0(rep(" ", 50 - n.bar), collapse = ""), "] - ",
      perc,
      "% - ETA - ",
      round(as.numeric((1 - raw.perc) * time.past / raw.perc, units = "mins"), 2),
      " minutes                 "
    ))

    structure <- do.call(function(...){set.block.value(ref.strucuture,...)},search.data[i, ])

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
    args <- outcomes
    args$structure <- structure
    args$smooth <- FALSE
    args$p.monit <- p.monit
    args$c.monit <- c.monit

    fitted.model <- do.call(fit_model, args)

    T_len <- fitted.model$t
    if (is.na(metric.cutoff)) {
      metric.cutoff <- floor(T_len / 10)
    }

    r <- length(fitted.model$outcomes)
    metric <- tolower(metric)
    predictions <- eval_past(fitted.model, eval_t = (metric.cutoff+1):T_len, lag = lag, pred.cred = pred.cred, eval.pred = TRUE)

    search.data$log.like[i] <- sum(predictions$log.like)
    search.data$mae[i] <- mean(predictions$mae)
    search.data$rae[i] <- mean(predictions$rae)
    search.data$mse[i] <- mean(predictions$mse)
    search.data$mase[i] <- mean(predictions$mase)
    search.data$interval.score[i] <- mean(predictions$interval.score)

    cur.metric <- ifelse(metric == "log.like", -search.data$log.like[i], search.data[[metric]][i])
    if (cur.metric < best.metric) {
      best.structure <- structure
      best.model <- fitted.model
      best.metric <- cur.metric
    }
  }
  cat("\n")
  search.data <- search.data[order(-search.data[[metric]],decreasing=(metric != "log.like")), ]
  if (smooth) {
    best.model <- smoothing(best.model)
  }

  out.vals <- list(
    search.data = search.data,
    structure = best.structure,
    model = best.model
  )
  class(out.vals) <- "searched_dlm"

  return(out.vals)
}

#' Auxiliary function for evaluating the posterior density of a DLM
#'
#' Evaluates the density for a set of parameters theta in a DLM. The structure of the DLM is taken to be that of the fitted_dlm object passed as input.
#'
#' @param theta Matrix: A matrix representing the set of parameter for which to evaluate the density. Its size should be n x t, where n is the number of latent states and t is the length of the time series;
#' @param model fitted_dlm: A fitted_dlm object.

#' @return A scalar representing the log density evaluated at theta.
#' @export
#'
#' @keywords internal
#' @examples
#'
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 0.95)
#' season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)
#'
#' outcome <- Poisson(lambda = "rate", data = data)
#'
#' fitted.data <- fit_model(level, season,
#'   AirPassengers = outcome
#' )
#' eval_dlm_post(fitted.data)
eval_dlm_post <- function(theta, model) {
  t <- model$t
  n <- model$n
  k <- model$k

  G <- model$G
  G.labs <- model$G.labs

  mts <- model$mt
  Cts <- model$Ct
  mt <- model$mt
  Ct <- model$Ct
  at <- model$at
  Rt <- model$Rt

  Ct.placeholder=Cts[,,1]*0

  log.post <- dmvnorm(theta[, t], model$mt[,t], Cts[, , t] |> matrix(n, n))

  if (t > 1) {
    for (i in (t - 1):1) {
      mt.step <- mt[, i, drop = FALSE]
      Ct.step <- Ct[, , i]
      Rt.step <- Rt[, , i+1]
      # G.step <- calc_current_G(mt.step, Ct.step, G[, , i + 1], G.labs)$G
      # G.step <- G[, , i + 1]
      G.step <- calc_current_G(theta[, i], Ct.placeholder, G[, , i + 1], G.labs)$G

      simple.Rt.inv <- Ct.step %*% transpose(G.step) %*% ginv(Rt.step)

      mts[, i] <- mt.step + simple.Rt.inv %*% (theta[, i + 1] - at[,i+1])
      Cts[, , i] <- Ct.step - simple.Rt.inv %*% Rt.step %*% transpose(simple.Rt.inv)
      # Cts[, , i] <- Ct.step - Ct.step %*% transpose(G.step) %*% ginv(Rt.step) %*% G.step %*% Ct.step
      log.post <- log.post + dmvnorm(theta[, i], mts[, i], Cts[, , i])
    }
  }
  return(log.post)
}

#' Auxiliary function for evaluating the prior density of a DLM
#'
#' Evaluates the prior density for a set of parameters theta in a DLM. The structure of the DLM is taken to be that of the fitted_dlm object passed as input.
#'
#' @param theta Matrix: A matrix representing the set of parameter for which to evaluate the density. Its size should be n x t, where n is the number of latent states and t is the length of the time series;
#' @param model fitted_dlm: A fitted_dlm object.
#'
#' @return A scalar representing the log density evaluated at theta.
#' @export
#'
#' @keywords internal
#' @examples
#'
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 0.95)
#' season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)
#'
#' outcome <- Poisson(lambda = "rate", data = data)
#'
#' fitted.data <- fit_model(level, season,
#'   AirPassengers = outcome
#' )
#' eval_dlm_prior(fitted.data)
eval_dlm_prior <- function(theta, model) {
  t <- model$t

  a1 <- model$a1
  R1 <- model$R1
  G <- model$G
  G.labs <- model$G.labs
  h <- model$h
  W <- model$W

  R1_placeholder <- R1 * 0
  D_placeholder <- R1**0

  log.prior <- dmvnorm(theta[, t], a1, R1)

  if (t > 1) {
    for (i in 2:t) {
      next.step <- one_step_evolve(theta[, i - 1], R1_placeholder, G[, , i], G.labs, D_placeholder, h[, i, drop = FALSE], W[, , i])

      log.prior <- log.prior + dmvnorm(theta[, i], next.step$at, next.step$Rt)
    }
  }
  return(log.prior)
}

#' Auxiliary function for evaluating the prior density of a DLM
#'
#' Evaluates the prior density for a set of parameters theta in a DLM. The structure of the DLM is taken to be that of the fitted_dlm object passed as input.
#'
#' @param theta Matrix: A matrix representing the set of parameter for which to evaluate the density. Its size should be n x t, where n is the number of latent states and t is the length of the time series;
#' @param model fitted_dlm: A fitted_dlm object.
#'
#' @return A scalar representing the log likelihood evaluated at theta.
#' @export
#'
#' @keywords internal
#' @examples
#'
#' data <- c(AirPassengers)
#'
#' level <- polynomial_block(rate = 1, order = 2, D = 0.95)
#' season <- harmonic_block(rate = 1, order = 2, period = 12, D = 0.975)
#'
#' outcome <- Poisson(lambda = "rate", data = data)
#'
#' fitted.data <- fit_model(level, season,
#'   AirPassengers = outcome
#' )
#' eval_dlm_log_like(fitted.data)
eval_dlm_log_like <- function(theta, model) {
  t <- model$t
  n <- model$n
  k <- model$k
  Ct.placeholder <- matrix(0, n, n)
  FF <- model$FF
  FF.labs <- model$FF.labs
  pred.names <- model$pred.names

  log.like <- 0
  for (i in seq_len(t)) {
    lin.pred <- calc_lin_pred(theta[, i, drop = FALSE], Ct.placeholder, FF[, , i] |> matrix(n, k), FF.labs, pred.names)
    ft <- lin.pred$ft
    Qt <- lin.pred$Qt
    for (outcome in model$outcomes) {
      offset.step <- outcome$offset[i, ]
      na.flag <- any(is.null(offset.step) | any(offset.step == 0) | any(is.na(offset.step)) | any(is.na(outcome$data[i, ])))
      pred.index <- match(outcome$pred.names, pred.names)
      ft.canom <- ft[pred.index, , drop = FALSE]
      Qt.canom <- Qt[pred.index, pred.index, drop = FALSE]
      if (outcome$convert.canom.flag) {
        ft.canom <- outcome$convert.mat.canom %*% ft.canom
        Qt.canom <- outcome$convert.mat.canom %*% Qt.canom %*% transpose(outcome$convert.mat.canom)
      }

      if (!na.flag) {
        offset.pred <- outcome$apply_offset(ft.canom, Qt.canom, offset.step)
        ft.canom <- offset.pred$ft
        Qt.canom <- offset.pred$Qt
      }
      param <- outcome$inv_link_function(ft.canom)
      log.like <- log.like + outcome$log_like_function(outcome$data[i, ], param)
    }
  }
  return(log.like)
}
