#' Fitting kDGLM models
#'
#' Fit a model given it's structure and the observed data. This function can be used for any supported family (see vignette).
#'
#' @param ... dlm_block object: The structural blocks of the model. All block must be completely defined.
#' @param outcomes List: The observed data. It should contain objects of the class dlm_distr.
#' @param pred_cred Numeric: A number between 0 and 1 (not included) indicating the credibility interval for predictions. If not within the valid range of values, predictions are not made.
#' @param smooth Bool: A flag indicating if the smoothed distribution for the latent variables should be calculated.
#' @param p_monit numeric (optional): The prior probability of changes in the latent space variables that are not part of it's dynamic.
#' @param c_monit numeric (optional, if p_monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
#'
#' @return The fitted model (a fitted_dlm object). Contains the values of the estimated parameter and some extra info regarding the quality of the fit.
#' @export
#'
#' @examples
#' library(kDGLM)
#'
#' # Poisson case
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
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
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
#' show_fit(fitted_data, smooth = TRUE)$plot
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
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Unknown variance
#' outcome <- Normal(mu = "mu", Sigma = "sigma2", outcome = data)
#'
#' fitted_data <- fit_model(level, variance, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' ##################################################################
#'
#' # Gamma case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' data <- matrix(rgamma(T, phi, phi / (20 * (sin(w * 1:T / T) + 2))), T, 1)
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
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # DO NOT RUN
#' # # Unknown shape
#' # outcome <- Gamma(phi = "phi", mu = "mu", outcome = data)
#' #
#' # fitted_data <- fit_model(level, season, scale, outcomes = outcome)
#' # summary(fitted_data)
#' #
#' # show_fit(fitted_data, smooth = TRUE)$plot
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

  structure <- block_merge(...)
  if (structure$status == "undefined") {
    stop("Error: One or more hiper parameter are undefined. Did you meant to use the search_model function?")
  }

  if (structure$t == 1) {
    structure$t <- t
    structure$G <- array(structure$G, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$G))
    structure$D <- array(structure$D, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$D))
    structure$W <- array(structure$W, c(structure$n, structure$n, structure$t), dimnames = dimnames(structure$W))
    structure$FF <- array(structure$FF, c(structure$n, structure$k, structure$t), dimnames = dimnames(structure$FF))
  }
  structure$G[, , 1] <- diag(structure$n)
  structure$D[, , 1] <- 1
  structure$W[, , 1] <- 0
  if (t != structure$t) {
    stop(paste0("Error: outcome does not have the same time length as structure: got ", t, " from outcome, expected ", structure$t))
  }

  var_names_out <- c()
  for (outcome in outcomes) {
    var_names_out <- c(var_names_out, outcome$var_names)
  }
  unique(var_names_out)
  if (any(!(var_names_out %in% structure$var_names))) {
    stop("Error: One or more linear predictor in outcomes are not present in the model structure.")
  }
  if (any(!(structure$var_names %in% var_names_out))) {
    warning("One or more linear predictor in the model structure are not used in the outcomes.")
  }

  model <- analytic_filter(
    outcomes = outcomes,
    m0 = structure$m0,
    C0 = structure$C0,
    FF = structure$FF,
    G = structure$G,
    G_labs = structure$G_labs,
    D = structure$D,
    W = structure$W,
    p_monit = p_monit,
    c_monit = c_monit
  )
  if (smooth) {
    smoothed <- generic_smoother(model$mt, model$Ct, model$at, model$Rt, model$G, model$G_labs)
    model$mts <- smoothed$mts
    model$Cts <- smoothed$Cts
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
  model$m0 <- structure$m0
  model$C0 <- structure$C0
  model$names <- structure$names
  model$smooth <- smooth
  model$pred_cred <- pred_cred
  model$t <- t
  class(model) <- "fitted_dlm"

  return(model)
}

#' forecast
#'
#' Makes predictions for t times ahead using the fitted model.
#'
#' @param model <undefined class> or list: The fitted model to be use for predictions.
#' @param t Numeric: Time window for prediction.
#' @param outcome List (optional): A named list containing the observed values in the prediction window. Note that the names in the list must be the same as the names passed during the fitting process.
#' @param offset Matrix or scalar: offset for predictions. Should have dimensions k x t, where k is the number of outcomes of the model. If offset is not specified, the last value observed by the model will be used.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x k x t, where n is the number of latent variables, k is the number of outcomes in the model. If not specified, the last value given to the model will be used.
#' @param G Array: A 3D-array containing the evolution matrix for each time. It's dimension should be n x n x t, where n is the number of latent variables. If not specified, the last value given to the model will be used.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x t, where n is the number of latent variables and T is the time series length. If not specified, the last value given to the model will be used in the first step, and 1 will be use thereafter.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be n x n x t, where n is the number of latent variables and T is the time series length. If not specified, 0 will be used.
#' @param plot Bool: A flag indicating if a plot should be produced.
#' @param pred_cred Numeric: The credibility level for the I.C. intervals.
#'
#' @return A list containing:
#' \itemize{
#'    \item pred Matrix: A matrix with the predictive mean at each time. Dimensions are k x t, where k is the number of outcomes.
#'    \item var.pred Array: A 3D-array with the predictive covariance matrix at each time. Dimensions are k x k x t, where k is the number of outcomes.
#'    \item icl.pred Matrix: A matrix with the lower bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#'    \item icu.pred Matrix: A matrix with the upper bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#'    \item at Matrix: A matrix with the values for the linear predictor at each time. Dimensions are k x t, where k is the number of outcomes.
#'    \item Rt Array: A 3D-array with the covariance of the linear predictor matrix at each time. Dimensions are k x k x t, where k is the number of outcomes.
#'    \item plot (if so chosen): A plotly or ggplot object.
#' }
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom Rfast transpose data.frame.to_matrix
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
#' forecast(fitted_data, 20)$plot
#'
#' @family {auxiliary functions for fitted_dlm objects}
forecast <- function(model, t = 1, outcome = NULL, offset = NULL, FF = NULL, G = NULL, D = NULL, W = NULL, plot = ifelse(requireNamespace("plotly", quietly = TRUE), "plotly", TRUE), pred_cred = 0.95) {
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  k <- dim(model$FF)[2]
  r <- sum(sapply(model$outcomes, function(x) {
    x$r
  }))
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
  if (length(dim(W)) > 3) {
    stop(paste0("Error: W should have at most 3 dimensions. Got ", length(dim(W)), "."))
  }
  if (length(dim(offset)) > 2) {
    stop(paste0("Error: D should have at most 2 dimensions. Got ", length(dim(offset)), "."))
  }

  G_labs <- model$G_labs
  if (is.null(G)) {
    G <- array(model$G[, , t_last], c(n, n, t))
  }
  if (is.null(FF)) {
    FF <- array(model$FF[, , t_last], c(n, k, t))
  }
  if (is.null(D)) {
    D <- array(model$D[, , t_last], c(n, n, t))
    D[, , -1] <- 1
  } else {
    if (all(dim(D) == 1)) {
      D <- array(D, c(n, n, t))
    } else {
      if (length(dim(D)) == 2 | (length(dim(D)) == 3 & dim(D)[3] == 1)) {
        D <- array(D, c(n, n, t))
        D[, , -1] <- 1
      }
    }
  }
  if (is.null(W)) {
    W <- array(model$W[, , t_last], c(n, n, t))
    W[, , -1] <- 0
  } else {
    if (all(dim(W) == 1)) {
      W <- array(diag(n) * W, c(n, n, t))
    } else {
      if (length(dim(W)) == 2 | (length(dim(W)) == 3 & dim(W)[3] == 1)) {
        W <- array(W, c(n, n, t))
        W[, , -1] <- 0
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
  if (any(dim(W) != c(n, n, t))) {
    stop(paste0("Error: W has wrong dimesions. Expected ", paste(n, n, t, sep = "x"), ". Got ", paste(dim(W), colapse = "x")))
  }
  #####

  m0 <- model$mt[, t_last]
  C0 <- model$Ct[, , t_last]

  m1 <- matrix(NA, n, t)
  C1 <- array(NA, c(n, n, t))
  f1 <- matrix(NA, k, t)
  Q1 <- array(NA, c(k, k, t))

  D <- ifelse(D == 0, 1, D)
  D_inv <- 1 / D

  last_m <- m0
  last_C <- C0

  outcome_forecast <- list()
  for (outcome_name in names(model$outcomes)) {
    r_i <- model$outcomes[[outcome_name]]$r
    outcome_forecast[[outcome_name]]$conj_param <- matrix(NA, t, length(model$outcomes[[outcome_name]]$conj_prior_param)) %>% as.data.frame()
    names(outcome_forecast[[outcome_name]]$conj_param) <- names(model$outcomes[[outcome_name]]$conj_prior_param)

    outcome_forecast[[outcome_name]]$ft <- matrix(NA, model$outcomes[[outcome_name]]$k, t)
    outcome_forecast[[outcome_name]]$Qt <- array(NA, c(model$outcomes[[outcome_name]]$k, model$outcomes[[outcome_name]]$k, t))

    if (!is.null(outcome)) {
      outcome_forecast[[outcome_name]]$outcome <- outcome[[outcome_name]] %>% matrix(t, r_i)
    } else {
      outcome_forecast[[outcome_name]]$outcome <- model$outcomes[[outcome_name]]$outcome[t_last, ] %>% matrix(t, r_i)
    }
    if (!is.null(offset)) {
      outcome_forecast[[outcome_name]]$offset <- offset[[outcome_name]] %>% matrix(t, r_i)
    } else {
      outcome_forecast[[outcome_name]]$offset <- model$outcomes[[outcome_name]]$offset[t_last, ] %>% matrix(t, r_i)
    }
  }


  for (t_i in c(1:t)) {
    next_step <- one_step_evolve(last_m, last_C, G[, , t_i] %>% matrix(n, n), G_labs, D_inv[, , t_i]**0, W[, , t_i] + C0 * (D_inv[, , t_i] - 1))
    last_m <- next_step$at
    last_C <- next_step$Rt

    m1[, t_i] <- last_m
    C1[, , t_i] <- last_C

    lin_pred <- calc_lin_pred(last_m, last_C, FF[, , t_i] %>% matrix(n, k))
    f1[, t_i] <- lin_pred$ft
    Q1[, , t_i] <- lin_pred$Qt
    for (outcome_name in names(model$outcomes)) {
      model_i <- model$outcomes[[outcome_name]]
      pred_index <- match(model_i$var_names, model$var_names)
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

  data_list <- list("obs" = output, "pred" = pred %>% data.frame.to_matrix() %>% t(), "icl.pred" = icl.pred %>% as.matrix() %>% t(), "icu.pred" = icu.pred %>% as.matrix() %>% t())
  data_name <- c("Observation", "Prediction", "C.I.lower", "C.I.upper")

  data_raw <- lapply(1:4, function(i) {
    data <- cbind(as.character(1:t + t_last) %>% as.data.frame(), data_list[[i]]) %>%
      as.data.frame() %>%
      pivot_longer(1:r + 1)

    data$name <- as.factor(data$name)
    names(data) <- c("Time", "Serie", data_name[i])
    levels(data$Serie) <- out_names
    data
  })
  data <- do.call(function(...) {
    data_list <- list(...)
    data <- data_list[[1]]
    for (i in 2:4) {
      data <- data %>% left_join(data_list[[i]], by = c("Time", "Serie"))
    }
    data
  }, data_raw)
  data$Time <- as.numeric(data$Time)
  data$Serie <- as.factor(data$Serie)

  data_past <- eval_past(model, smooth = TRUE, pred_cred = pred_cred)
  plot_data <- rbind(
    cbind(data_past, type = "Past"),
    cbind(data, type = "Forecast")
  )

  return_list <- list(data = plot_data, forecast = data, outcomes = outcome_forecast, mt = m1, Ct = C1, ft = f1, Qt = Q1)

  if (plot != FALSE) {
    obs_na_rm <- data_past$Observation[!is.na(data_past$Observation)]
    max_value <- calcula_max(obs_na_rm - min(obs_na_rm))[[3]] + min(obs_na_rm)
    min_value <- -calcula_max(-(obs_na_rm - max(obs_na_rm)))[[3]] + max(obs_na_rm)
    plot_data$shape_point <- paste0("Obs_", plot_data$type)
    plot_data$group_ribbon <- paste0(plot_data$Serie, plot_data$type)


    plt.obj <- ggplot(plot_data, aes_string(x = "Time")) +
      geom_line(aes_string(y = "Prediction", linetype = "type", color = "Serie")) +
      geom_ribbon(aes_string(ymin = "C.I.lower", ymax = "C.I.upper", fill = "Serie", group = "group_ribbon"), alpha = 0.25, color = NA) +
      geom_point(aes_string(y = "Observation", shape = "shape_point", color = "Serie")) +
      scale_fill_hue("", na.value = NA) +
      scale_color_hue("", na.value = NA) +
      scale_linetype_manual("", values = c("dashed", "solid")) +
      scale_shape_manual("", values = c(17, 16)) +
      scale_x_continuous("Time") +
      scale_y_continuous("$y_t$") +
      theme_bw() +
      coord_cartesian(ylim = c(min_value, max_value))
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

  return(return_list)
}

#' eval_past
#'
#' Evaluates the predictive values for the observed values used to fit the model. Predictions can be made with smoothed values or with filtered values with a time offset.
#'
#' @param model <undefined class> or list: The fitted model to be use for evaluation.
#' @param smooth bool: The flag indicating if smoothed values should be used. If TRUE, h will not be used.
#' @param h positive integer: The relative offset for forecast. Values for time t will be calculated based on the filtered values of time t-h Will be ignored if smooth is TRUE.
#' @param pred_cred Numeric: The credibility level for the I.C. intervals.
#'
#' @return A list containing:
#' \itemize{
#'    \item pred Matrix: A matrix with the predictive mean at each time. Dimensions are k x t, where k is the number of outcomes.
#'    \item var.pred Array: A 3D-array with the predictive covariance matrix at each time. Dimensions are k x k x t, where k is the number of outcomes.
#'    \item icl.pred Matrix: A matrix with the lower bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#'    \item icu.pred Matrix: A matrix with the upper bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
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
#' past <- eval_past(fitted_data, smooth = TRUE)
#'
#' @family {auxiliary functions for fitted_dlm objects}
eval_past <- function(model, smooth = FALSE, h = 0, pred_cred = 0.95) {
  if (smooth & h > 0) {
    h <- 0
    warning("h is only used if smooth is set to TRUE.")
  }
  if (h < 0 | round(h) != h) {
    stop(paste0("ERROR: h should be a positive integer. Got ", h, "."))
  }
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  r <- dim(model$FF)[2]
  k <- sum(sapply(model$outcomes, function(x) {
    x$r
  }))

  FF <- model$FF
  G <- array(diag(n), c(n, n, t_last))
  pred <- matrix(NA, k, t_last)
  var.pred <- array(NA, c(k, k, t_last))
  icl.pred <- matrix(NA, k, t_last)
  icu.pred <- matrix(NA, k, t_last)
  log.like <- rep(NA, t_last)

  if (smooth) {
    ref_mt <- model$mts
    ref_Ct <- model$Cts
    D <- model$D**0
    D_inv <- 1 / D
    W <- model$W * 0
  } else {
    ref_mt <- model$mt
    ref_Ct <- model$Ct
    D <- model$D
    D_inv <- 1 / D
    W <- model$W
    if (h > 0) {
      # for (i in c(1:h)) {
      #   G <- G %*% model$G
      # }
      G <- model$G
    }
  }
  G_labs <- model$G_labs

  for (i in c((1 + h):t_last)) {
    mt <- if (i <= h) {
      model$m0
    } else {
      ref_mt[, (i - h):(i - h)]
    }
    Ct <- if (i <= h) {
      model$C0
    } else {
      ref_Ct[, , i - h]
    }
    # model$filter(model$outcome[i, ], mt, Ct, FF[, , i] %>% matrix(n, r), G, D_inv[, , i], W[, , i], model$offset[i, ], parms = model$parms)
    next_step <- list("at" = mt, "Rt" = Ct)
    if (h >= 1) {
      for (t in c(1:h)) {
        next_step <- one_step_evolve(next_step$at, next_step$Rt, G[, , i - t + 1], G_labs, D_inv[, , i - t + 1]**(t == 1), W[, , i - t + 1])
      }
    }
    # next_step <- one_step_evolve(mt, Ct, FF[, , i] %>% matrix(n, r), G, G_labs, D_inv[, , i], W[, , i])
    lin_pred <- calc_lin_pred(next_step$at %>% matrix(n, 1), next_step$Rt, FF[, , i] %>% matrix(n, r))

    r_acum <- 0
    for (outcome in model$outcomes) {
      r_cur <- outcome$r

      pred_index <- match(outcome$var_names, model$var_names)
      k_i <- length(pred_index)

      cur_step <- outcome$apply_offset(lin_pred$ft[pred_index, , drop = FALSE], lin_pred$Qt[pred_index, pred_index, drop = FALSE], outcome$offset[i, ])

      if (outcome$convert_canom_flag) {
        ft_canom <- outcome$convert_mat_canom %*% cur_step$ft
        Qt_canom <- outcome$convert_mat_canom %*% cur_step$Qt %*% transpose(outcome$convert_mat_canom)
      } else {
        ft_canom <- cur_step$ft
        Qt_canom <- cur_step$Qt
      }


      conj_distr <- outcome$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)
      prediction <- outcome$calc_pred(conj_distr, outcome$outcome[i, , drop = FALSE], pred_cred, parms = outcome$parms)

      pred[(r_acum + 1):(r_acum + r_cur), i] <- prediction$pred
      var.pred[(r_acum + 1):(r_acum + r_cur), (r_acum + 1):(r_acum + r_cur), i] <- prediction$var.pred
      icl.pred[(r_acum + 1):(r_acum + r_cur), i] <- prediction$icl.pred
      icu.pred[(r_acum + 1):(r_acum + r_cur), i] <- prediction$icu.pred
      log.like[i] <- log.like[i] + sum(prediction$log.like)
      r_acum <- r_acum + r_cur
    }
  }

  r_acum <- 0
  out_names <- rep(NA, k)
  output <- matrix(NA, t_last, k)
  for (outcome_name in names(model$outcomes)) {
    r_cur <- model$outcomes[[outcome_name]]$r
    if (r_cur > 1) {
      out_names[(r_acum + 1):(r_acum + r_cur)] <- paste0(outcome_name, "_", 1:r_cur)
    } else {
      out_names[(r_acum + 1):(r_acum + r_cur)] <- outcome_name
    }
    output[, (r_acum + 1):(r_acum + r_cur)] <- model$outcomes[[outcome_name]]$outcome

    r_acum <- r_acum + r_cur
  }

  data_list <- list("obs" = output, "pred" = pred %>% data.frame.to_matrix() %>% t(), "icl.pred" = icl.pred %>% as.matrix() %>% t(), "icu.pred" = icu.pred %>% as.matrix() %>% t())
  data_name <- c("Observation", "Prediction", "C.I.lower", "C.I.upper")

  data_raw <- lapply(1:4, function(i) {
    data <- cbind(as.character(1:t_last) %>% as.data.frame(), data_list[[i]]) %>%
      as.data.frame() %>%
      pivot_longer(1:k + 1)

    data$name <- as.factor(data$name)
    names(data) <- c("Time", "Serie", data_name[i])
    levels(data$Serie) <- out_names
    data
  })
  data <- do.call(function(...) {
    data_list <- list(...)
    data <- data_list[[1]]
    for (i in 2:4) {
      data <- data %>% left_join(data_list[[i]], by = c("Time", "Serie"))
    }
    data
  }, data_raw)
  data$Time <- as.numeric(data$Time)
  data$Serie <- as.factor(data$Serie)
  return(data)
}

#' dlm_sampling
#'
#' This is function draws samples from the latent states using the FFBS algorithm. See \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 15, for details.
#'
#' @param model fitted_dlm: A fitted model from which to sample.
#' @param sample_size integer: The number of samples to draw.
#' @param filtered_distr bool: A flag indicating if the samples should be dawned from the filtered (if TRUE) or smoothed distribution. The default is FALSE.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Array: An array containing the samples of latent states. Dimensions are n x T x sample_size, where n is the number of latent variable in the model and T is the number of observed values.
#'    \item ft Array: An array containing the samples of linear predictors. Dimensions are m x T x sample_size, where m is the number of linear predictors in the model and T is the number of observed values.
#'    \item param Array: An array containing the samples of the parameters of the observational model. Dimensions are k x T x sample_size, where k is the number of parameters in the observational model and T is the number of observed values.
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
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' sample <- dlm_sampling(fitted_data, 2000)
#'
#' @family {auxiliary functions for fitted_dlm objects}
dlm_sampling <- function(model, sample_size, filtered_distr = FALSE) {
  G <- model$G
  G_labs <- model$G_labs
  at <- model$at
  mt <- model$mt
  mts <- model$mts
  FF <- model$FF
  T_len <- dim(mt)[2]
  n <- dim(mt)[1]
  k <- dim(FF)[2]
  outcomes <- list()
  for (outcome_name in names(model$outcomes)) {
    outcomes[[outcome_name]] <- list(
      inv_link = model$outcomes[[outcome_name]]$inv_link_function,
      apply_offset = model$outcomes[[outcome_name]]$apply_offset,
      offset = model$outcomes[[outcome_name]]$offset,
      l = model$outcomes[[outcome_name]]$l,
      convert_mat_canom = model$outcomes[[outcome_name]]$convert_mat_canom,
      convert_canom_flag = model$outcomes[[outcome_name]]$convert_canom_flag
    )
  }

  mt_sample <- rnorm(n * T_len * sample_size) %>% array(c(n, T_len, sample_size))
  ft_sample <- array(NA, c(k, T_len, sample_size))

  Ct_chol <- var_decomp(model$Ct[, , T_len])
  # mt_sample_i <- transpose(Ct_chol) %*% mt_sample[, T_len, ] + mt[, T_len]
  mt_sample_i <- crossprod(Ct_chol, mt_sample[, T_len, ]) + mt[, T_len]
  FF_step <- FF[, , T_len]
  if (any(is.na(FF_step))) {
    ft_sample_i <- sapply(1:sample_size,
      function(j) {
        calc_lin_pred(mt_sample_i[, j], Ct_chol * 0, FF_step)$ft
      },
      simplify = "matrix"
    )
  } else {
    ft_sample_i <- crossprod(FF_step, mt_sample_i)
  }

  var_names <- colnames(FF)
  for (outcome_name in names(outcomes)) {
    pred_index <- match(model$outcome[[outcome_name]]$var_names, var_names)
    k_i <- length(pred_index)

    outcomes[[outcome_name]]$k_i <- k_i
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

      G_ref <- calc_current_G(mt_now, G[, , t + 1], G_labs)$G
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
          calc_lin_pred(mt_sample_i[, j], Ct_chol * 0, FF_step)$ft
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
