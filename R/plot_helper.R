#' Visualizing a fitted kDGLM model
#'
#' Calculate the predictive mean and some quantile for the observed data and show a plot.
#'
#' @param model fitted_dlm: A fitted DGLM.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param lag Integer: A integer with the number of steps ahead should be used for prediction. If lag<0, the smoothed distribution is used and, if lag==0, the filtered distribuition is used.
#' @param plotly Bool: A flag indicating if plotly should be used for creating plots.
#'
#' @return A list containing:
#' \itemize{
#'    \item data tibble object: A data frame containing the observations, predictions and credibility intervals at each time.
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#' }
#' @export
#' @importFrom grDevices rainbow
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
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
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
show_fit <- function(model, pred_cred = 0.95, lag = 1, plotly = requireNamespace("plotly", quietly = TRUE)) {
  t_last <- dim(model$mt)[2]
  eval <- eval_past(model, lag = lag, pred_cred = pred_cred, eval_pred = TRUE)$data

  obs_na_rm <- eval$Observation[!is.na(eval$Observation)]
  max_value <- evaluate_max(obs_na_rm - min(obs_na_rm))[[3]] + min(obs_na_rm)
  min_value <- -evaluate_max(-(obs_na_rm - max(obs_na_rm)))[[3]] + max(obs_na_rm)

  n_colors <- length(unique(eval$Serie))
  colors <- rainbow(n_colors, s = 0.6)
  series_names <- unique(eval$Serie)
  names(colors) <- series_names
  colors[["Detected changes"]] <- "black"
  linetypes <- c(
    "Detected changes" = "dashed",
    "Observation" = NA,
    "Fitted values" = "solid"
  )
  shapes <- c(
    "Detected changes" = NA,
    "Observation" = 16,
    "Fitted values" = NA
  )

  title <- if (lag < 0) {
    "Smoothed predictions"
  } else if (lag == 0) {
    "Filtered predictions"
  } else if (lag == 1) {
    "One-step-ahead predictions"
  } else {
    paste0(lag, "-steps-ahead predictions")
  }

  plt <- ggplot() +
    geom_ribbon(data = eval, aes_string(x = "Time", fill = "Serie", ymin = "C.I.lower", ymax = "C.I.upper"), alpha = 0.25) +
    geom_line(data = eval, aes_string(x = "Time", color = "Serie", y = "Prediction", linetype = "'Fitted values'")) +
    geom_point(data = eval, aes_string(x = "Time", color = "Serie", y = "Observation", shape = '"Observation"'), alpha = 0.5) +
    scale_linetype_manual("", values = linetypes) +
    scale_shape_manual("", values = shapes) +
    scale_fill_manual("", na.value = NA, values = colors) +
    scale_color_manual("", na.value = NA, values = colors) +
    scale_y_continuous(name = "$y_t$") +
    scale_x_continuous("Time") +
    ggtitle(title) +
    theme_bw() +
    coord_cartesian(ylim = c(min_value, max_value))
  for (name_i in names(model$outcomes)) {
    outcome_i <- model$outcomes[[name_i]]
    if (any(outcome_i$alt.flags == 1)) {
      plt <- plt +
        geom_vline(
          data = data.frame(xintercept = (1:t_last)[outcome_i$alt.flags == 1], linetype = "Detected changes", serie = name_i),
          aes_string(xintercept = "xintercept", linetype = "linetype")
        )
    }
  }
  if (plotly) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("The plotly package is required for plotly plots.")
    } else {
      plt <- plotly::ggplotly(plt)

      for (i in (1:n_colors) - 1) {
        plt$x$data[[i + 1]]$legendgroup <-
          plt$x$data[[i + 1 + n_colors]]$legendgroup <-
          plt$x$data[[i + 1]]$name <-
          plt$x$data[[i + 1 + n_colors]]$name <- paste0(series_names[i + 1], ": fitted values")

        plt$x$data[[i + 1]]$showlegend <- FALSE

        plt$x$data[[i + 1 + 2 * n_colors]]$legendgroup <-
          plt$x$data[[i + 1 + 2 * n_colors]]$name <- paste0(series_names[i + 1], ": observations")
      }
      n <- length(plt$x$data)
      if (n %% 3 == 1) {
        plt$x$data[[n]]$legendgroup <-
          plt$x$data[[n]]$name <- "Detected changes"
      }
    }
  }
  return(list("data" = eval, "plot" = plt))
}

#' Visualizing latent states in a fitted kDGLM model
#'
#' @param model fitted_dlm: A fitted DGLM model.
#' @param var Character: The name of the variables to plot (same value passed while creating the structure). Any variable whose name partially match this variable will be plotted.
#' @param lag Integer: A integer with the number of steps ahead should be used for evaluating the latent variables. If lag<0, the smoothed distribution is used and, if lag==0, the filtered distribution is used.
#' @param cutoff Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param plotly Bool: A flag indicating if plotly should be used to create plots.
#'
#' @return A list containing:
#' \itemize{
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#'    \item data data.frame: A data frame containing the data used in the plot.
#' }
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom grDevices rainbow
#' @importFrom stats qnorm
#'
#' @examples
#'
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
#' plot_lat_var(fitted_data, smooth = TRUE)$plot
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
plot_lat_var <- function(model, var = "", lag = -1, cutoff = floor(model$t / 10), pred_cred = 0.95, plotly = requireNamespace("plotly", quietly = TRUE)) {
  if (!any(grepl(var, model$var_labels))) {
    stop(paste0("Error: Invalid variable selection. Got '", var, "', expected one of the following:\n", paste0(names(model$var_names), collapse = "\n")))
  }
  if (pred_cred >= 1 | pred_cred <= 0) {
    stop(paste0("Error: Invalid value for I.C. width. Must be between 0 and 1, got ", pred_cred))
  }

  indice <- (1:model$n)[grepl(var, model$var_labels)]
  var_names <- model$var_labels[grepl(var, model$var_labels)]

  size <- length(indice)
  t <- model$t

  param_distr <- eval_past(model, lag = lag, pred_cred = pred_cred, eval_pred = FALSE)

  m1 <- param_distr$mt[indice, (cutoff + 1):t, drop = FALSE] %>%
    t()
  std_mat <- param_distr$Ct[indice, indice, (cutoff + 1):t, drop = FALSE] %>%
    apply(3, diag) %>%
    sqrt() %>%
    matrix(size, t - cutoff) %>%
    t()

  lim_i <- m1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- m1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  m1 <- as.data.frame(m1)
  lim_i <- as.data.frame(lim_i)
  lim_s <- as.data.frame(lim_s)

  names(m1) <- var_names
  names(lim_i) <- var_names
  names(lim_s) <- var_names

  max_value <- evaluate_max(m1 - min(m1))[[3]] + min(m1)
  min_value <- -evaluate_max(-(m1 - max(m1)))[[3]] + max(m1)

  if (max_value - min_value < 1e-2) {
    center <- (max_value + min_value) / 2
    max_value <- center + 0.01
    max_value <- center - 0.01
  }


  m1$time <- c(1:dim(m1)[1]) + cutoff
  lim_i$time <- c(1:dim(lim_i)[1]) + cutoff
  lim_s$time <- c(1:dim(lim_s)[1]) + cutoff

  m1 <- m1 %>%
    pivot_longer(1:size) %>%
    rename("media" = "value")
  lim_i <- lim_i %>%
    pivot_longer(1:size) %>%
    rename("lim_i" = "value")
  lim_s <- lim_s %>%
    pivot_longer(1:size) %>%
    rename("lim_s" = "value")

  label <- paste0("\n(cred. ", pred_cred * 100 %>% round(), "%)")

  plot_data <- m1 %>%
    inner_join(lim_i, by = c("time", "name")) %>%
    inner_join(lim_s, by = c("time", "name"))

  plot_data$name <- factor(plot_data$name, levels = unique(plot_data$name))

  n_var <- length(levels(plot_data$name))
  color_list <- rainbow(n_var, s = 0.5)
  names(color_list) <- paste(levels(plot_data$name), label)

  fill_list <- rainbow(n_var, s = 0.5)
  names(fill_list) <- paste(levels(plot_data$name), label)
  plot_data$fill_name <- paste(plot_data$name, label)
  plot_data$color_name <- paste(plot_data$name, label)

  title <- if (lag < 0) {
    "Smoothed estimation of latent variables"
  } else if (lag == 0) {
    "Filtered estimation of latent variables"
  } else if (lag == 1) {
    "One-step-ahead prediction for latent variables"
  } else {
    paste0(lag, "-steps-ahead prediction for latent variables")
  }

  plt <- ggplot(plot_data, aes_string(x = "time", fill = "fill_name", color = "color_name")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous("Time") +
    scale_color_manual("", values = color_list, na.value = NA) +
    scale_fill_manual("", values = fill_list, na.value = NA) +
    labs(title = title) +
    scale_y_continuous("Parameter value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_ribbon(aes_string(ymin = "lim_i", ymax = "lim_s"), alpha = 0.25, color = NA) +
    geom_line(aes_string(y = "media")) +
    coord_cartesian(ylim = c(min_value, max_value))
  if (plotly) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("The plotly package is required for plotly plots.")
    } else {
      plt <- plotly::ggplotly(plt)

      for (i in (1:size) - 1) {
        plt$x$data[[i + 1 + 1]]$legendgroup <-
          plt$x$data[[i + 1 + size + 1]]$legendgroup <-
          plt$x$data[[i + 1 + 1]]$name <-
          plt$x$data[[i + 1 + size + 1]]$name <- plt$x$data[[i + 1 + size + 1]]$name

        plt$x$data[[i + 1 + 1]]$showlegend <- FALSE
      }
    }
  }
  return(list("plot" = plt, "data" = plot_data))
}

#' Visualizing linear predictors in a fitted kDGLM model
#'
#' @param model fitted_dlm: A fitted DGLM model.
#' @param pred Character: The name of the linear predictors to plot (same value passed while creating the structure). Any predictors whose name partially match this variable will be plotted.
#' @param lag Integer: A integer with the number of steps ahead should be used for evaluating the linear predictors. If lag<0, the smoothed distribution is used and, if lag==0, the filtered distribution is used.
#' @param cutoff Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param plotly Bool: A flag indicating if plotly should be used to create plots.
#'
#' @return A list containing:
#' \itemize{
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#'    \item data tibble: A data frame containing the data used in the plot.
#' }
#' @export
#' @importFrom grDevices rainbow
#' @importFrom stats qnorm
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#'
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
#' plot_lin_pred(fitted_data, smooth = TRUE)$plot
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
plot_lin_pred <- function(model, pred = "", lag = -1, cutoff = floor(model$t / 10), pred_cred = 0.95, plotly = requireNamespace("plotly", quietly = TRUE)) {
  if (!any(grepl(pred, model$pred_names))) {
    stop(paste0("Error: Invalid selected variable. Got ", pred, ", expected one of the following:\n", paste0(names(model$pred_names), collapse = "\n")))
  }
  if (pred_cred >= 1 | pred_cred <= 0) {
    stop(paste0("Error: Invalid value for I.C. width. Must be between 0 and 1, got ", pred_cred))
  }

  indice <- (1:model$k)[grepl(pred, model$pred_names)]
  var_names <- model$pred_names[grepl(pred, model$pred_names)]

  size <- length(indice)
  t <- model$t

  param_distr <- eval_past(model, lag = lag, pred_cred = pred_cred, eval_pred = FALSE)

  f1 <- param_distr$ft[indice, (cutoff + 1):t, drop = FALSE] %>%
    t()
  std_mat <- param_distr$Qt[indice, indice, (cutoff + 1):t, drop = FALSE] %>%
    apply(3, diag) %>%
    sqrt() %>%
    matrix(size, t - cutoff) %>%
    t()

  lim_i <- f1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- f1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  f1 <- as.data.frame(f1)
  lim_i <- as.data.frame(lim_i)
  lim_s <- as.data.frame(lim_s)

  names(f1) <- var_names
  names(lim_i) <- var_names
  names(lim_s) <- var_names

  max_value <- evaluate_max(f1 - min(f1))[[3]] + min(f1)
  min_value <- -evaluate_max(-(f1 - max(f1)))[[3]] + max(f1)

  if (max_value - min_value < 1e-2) {
    center <- (max_value + min_value) / 2
    max_value <- center + 0.01
    max_value <- center - 0.01
  }

  f1$time <- c(1:dim(f1)[1]) + cutoff
  lim_i$time <- c(1:dim(lim_i)[1]) + cutoff
  lim_s$time <- c(1:dim(lim_s)[1]) + cutoff

  f1 <- f1 %>%
    pivot_longer(1:size) %>%
    rename("media" = "value")
  lim_i <- lim_i %>%
    pivot_longer(1:size) %>%
    rename("lim_i" = "value")
  lim_s <- lim_s %>%
    pivot_longer(1:size) %>%
    rename("lim_s" = "value")

  label <- paste0("\n(cred. ", pred_cred * 100 %>% round(), "%)")

  plot_data <- f1 %>%
    inner_join(lim_i, by = c("time", "name")) %>%
    inner_join(lim_s, by = c("time", "name"))

  n_var <- length(unique(plot_data$name))
  color_list <- rainbow(n_var, s = 0.5)
  names(color_list) <- paste(unique(plot_data$name), label)

  fill_list <- rainbow(n_var, s = 0.5)
  names(fill_list) <- paste(unique(plot_data$name), label)
  plot_data$fill_name <- paste(plot_data$name, label)
  plot_data$color_name <- paste(plot_data$name, label)
  plot_data$IC_name <- paste(plot_data$name, label)

  title <- if (lag < 0) {
    "Smoothed estimation of linear predictors"
  } else if (lag == 0) {
    "Filtered estimation of linear predictors"
  } else if (lag == 1) {
    "One-step-ahead prediction for linear predictors"
  } else {
    paste0(lag, "-steps-ahead prediction for linear predictors")
  }

  plt <- ggplot(plot_data, aes_string(x = "time", fill = "fill_name", color = "color_name")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous("Time") +
    scale_color_manual("", values = color_list, na.value = NA) +
    scale_fill_manual("", values = fill_list, na.value = NA) +
    labs(title = title) +
    scale_y_continuous("predictor value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_ribbon(aes_string(ymin = "lim_i", ymax = "lim_s"), alpha = 0.25, color = NA) +
    geom_line(aes_string(y = "media")) +
    coord_cartesian(ylim = c(min_value, max_value))
  if (plotly) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("The plotly package is required for plotly plots.")
    } else {
      plt <- plotly::ggplotly(plt)

      for (i in (1:size) - 1) {
        plt$x$data[[i + 1 + 1]]$legendgroup <-
          plt$x$data[[i + 1 + size + 1]]$legendgroup <-
          plt$x$data[[i + 1 + 1]]$name <-
          plt$x$data[[i + 1 + size + 1]]$name <- plt$x$data[[i + 1 + size + 1]]$name

        plt$x$data[[i + 1 + 1]]$showlegend <- FALSE
      }
    }
  }
  return(list("plot" = plt, "data" = plot_data))
}
