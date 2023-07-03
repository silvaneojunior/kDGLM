#' Visualizing a fitted kDGLM model
#'
#' Calculate the predictive mean and some quantile for the observed data and show a plot.
#'
#' @param model fitted_dlm: A fitted DGLM.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param smooth Bool: A flag indicating if the smoothed should be used. If false, the filtered distribution will be used.
#' @param plotly Bool: A flag indicating if plotly should be used for creating plots.
#' @param h Integer: A integer with the amount of steps ahead should be used for prediction. Only used if smooth is false.
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
show_fit <- function(model, pred_cred = 0.95, smooth = model$smooth, plotly = requireNamespace("plotly", quietly = TRUE), h = 0) {
  t_last <- dim(model$mt)[2]
  eval <- eval_past(model, smooth = smooth, h = h, pred_cred = pred_cred)

  obs_na_rm <- eval$Observation[!is.na(eval$Observation)]
  max_value <- calcula_max(obs_na_rm - min(obs_na_rm))[[3]] + min(obs_na_rm)
  min_value <- -calcula_max(-(obs_na_rm - max(obs_na_rm)))[[3]] + max(obs_na_rm)

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
    theme_bw() +
    coord_cartesian(ylim = c(min_value, max_value))

  if (any(model$outcomes[[1]]$alt.flags == 1)) {
    plt <- plt +
      geom_vline(
        data = data.frame(xintercept = (1:t_last)[model$outcomes[[1]]$alt.flags == 1], linetype = "Detected changes"),
        aes_string(xintercept = "xintercept", linetype = "linetype", fill = "linetype", color = "linetype")
      )
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
    }
  }
  return(list("data" = eval, "plot" = plt))
}

#' Visualizing latent states in a fitted kDGLM model
#'
#' @param model <undefined class> or list: A fitted DGLM model.
#' @param var Character: The name of the variables to plot (same value passed while creating the structure). Any variable whose name partially match this variable will be plotted.
#' @param smooth Bool: A flag indicating if the smoothed distribution should be used. If false, the filtered distribution shall be used.
#' @param cut_off Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param plotly Bool: A flag indicating if plotly should be used to create plots.
#'
#' @return A list containing:
#' \itemize{
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#'    \item data tibble: A data frame containing the data used in the plot.
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
plot_lat_var <- function(model, var = "", smooth = model$smooth, cut_off = round(model$t/10), pred_cred = 0.95, plotly = requireNamespace("plotly", quietly = TRUE)) {
  if (!any(grepl(var, names(model$names)))) {
    stop(paste0("Error: Invalid selected variable. Got ", var, ", expected one of the following:\n", names(model$names)))
  }
  if (pred_cred >= 1 | pred_cred <= 0) {
    stop(paste0("Error: Invalid value for I.C. width. Must be between 0 and 1, got ", pred_cred))
  }

  indice <- c()
  names_var <- c()
  for (i in names(model$names)) {
    if (grepl(var, i)) {
      indice <- c(indice, model$names[[i]])
      count <- 0
      k_i=length(model$names[[i]])
      size=floor(log10(k_i))+1
      size_i=size-floor(log10(1:k_i))-1
      names_index=sapply(size_i,function(x){
        paste0(rep('0',x),collapse='')
      })
      for (index in model$names[[i]]) {
        count <- count + 1
        if (length(model$names[[i]]) > 1) {
          names_var <- c(names_var, paste0(i, ":\nlat. val. ",names_index[count], count))
        } else {
          names_var <- c(names_var, i)
        }
      }
    }
  }
  size <- length(indice)
  t <- dim(model$mts)[2]
  m1 <- if (smooth) {
    model$mts[indice, ]
  } else {
    model$mt[indice, ]
  }
  m1 <- m1 %>%
    matrix(size, t) %>%
    t()
  std_mat <- if (smooth) {
    model$Cts[indice, indice, ]
  } else {
    model$Ct[indice, indice, ]
  }
  if (size > 1) {
    std_mat <- std_mat %>% apply(3, diag)
  }
  std_mat <- std_mat %>%
    sqrt() %>%
    matrix(size, t) %>%
    t()

  lim_i <- m1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- m1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  m1 <- as.data.frame(m1)
  lim_i <- as.data.frame(lim_i)
  lim_s <- as.data.frame(lim_s)

  names(m1) <- names_var
  names(lim_i) <- names_var
  names(lim_s) <- names_var

  max_value <- calcula_max(m1 - min(m1))[[3]] + min(m1)
  min_value <- -calcula_max(-(m1 - max(m1)))[[3]] + max(m1)

  if (max_value - min_value < 1e-2) {
    center <- (max_value + min_value) / 2
    max_value <- center + 0.01
    max_value <- center - 0.01
  }

  m1$time <- c(1:dim(m1)[1])
  lim_i$time <- c(1:dim(lim_i)[1])
  lim_s$time <- c(1:dim(lim_s)[1])

  m1 <- m1[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename("media" = "value")
  lim_i <- lim_i[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename("lim_i" = "value")
  lim_s <- lim_s[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename("lim_s" = "value")

  label <- paste0("\n(", pred_cred * 100 %>% round(), "%)")

  plot_data <- m1 %>%
    inner_join(lim_i, by = c("time", "name")) %>%
    inner_join(lim_s, by = c("time", "name"))

  n_var <- length(unique(plot_data$name))
  color_list <- rainbow(n_var, s = 0.5)
  names(color_list) <- paste(unique(plot_data$name), label)

  fill_list <- rainbow(n_var, s = 0.5)
  names(fill_list) <- paste(unique(plot_data$name), label)
  plot_data$fill_name <- paste(plot_data$name, label)
  plot_data$color_name <- paste(plot_data$name, label)

  plt <- ggplot(plot_data, aes_string(x = "time", fill = "fill_name", color = "color_name")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous("Time") +
    scale_color_manual("", values = color_list, na.value = NA) +
    scale_fill_manual("", values = fill_list, na.value = NA) +
    labs(title = paste0(var, " (", ifelse(smooth, "smoothed", "only filtered"), ")")) +
    scale_y_continuous("Parameter value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_ribbon(aes_string(ymin = "lim_i", ymax = "lim_s"), alpha = 0.25) +
    geom_line(aes_string(y = "media")) +
    coord_cartesian(ylim = c(min_value, max_value))
  if (plotly) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("The plotly package is required for plotly plots.")
    } else {
      plt <- plotly::ggplotly(plt)

      # for (i in (1:size) - 1) {
      #   plt$x$data[[i + 1+1]]$legendgroup <-
      #     plt$x$data[[i + 1 + size+1]]$legendgroup <-
      #     plt$x$data[[i + 1+1]]$name <-
      #     plt$x$data[[i + 1 + size+1]]$name <- plt$x$data[[i + 1 + size+1]]$name
      #
      #   plt$x$data[[i + 1+1]]$showlegend <- FALSE
      # }
    }
  }
  return(list("plot" = plt, "data" = plot_data))
}

#' Visualizing linear predictors in a fitted kDGLM model
#'
#' @param model <undefined class> or list: A fitted DGLM model.
#' @param pred Character: The name of the linear predictors to plot (same value passed while creating the structure). Any predictors whose name partially match this variable will be plotted.
#' @param smooth Bool: A flag indicating if the smoothed distribution should be used. If false, the filtered distribution shall be used.
#' @param cut_off Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.
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
plot_lin_pred <- function(model, pred = "", smooth = model$smooth, cut_off = 10, pred_cred = 0.95, plotly = requireNamespace("plotly", quietly = TRUE)) {
  if (!any(grepl(pred, model$var_names))) {
    stop(paste0("Error: Invalid selected variable. Got ", pred, ", expected one of the following:\n", model$var_names))
  }
  if (pred_cred >= 1 | pred_cred <= 0) {
    stop(paste0("Error: Invalid value for I.C. width. Must be between 0 and 1, got ", pred_cred))
  }

  indice <- (1:length(model$var_names))[grepl(pred, model$var_names)]
  names_var <- model$var_names[grepl(pred, model$var_names)]

  size <- length(indice)
  t <- dim(model$mts)[2]
  m1 <- if (smooth) {
    model$mts
  } else {
    model$mt
  }
  C1 <- if (smooth) {
    model$Cts
  } else {
    model$Ct
  }
  f1 <- (sapply(1:t, function(t) {
    t(model$FF[, indice, t]) %*% m1[, t]
  }) %>% matrix(size, t))[, ]
  R1 <- (sapply(1:t, function(t) {
    sqrt(diag(t(model$FF[, indice, t]) %*% C1[, , t] %*% model$FF[, indice, t]))
  }) %>% matrix(size, t))
  f1 <- f1 %>%
    matrix(size, t) %>%
    t()
  std_mat <- R1 %>%
    sqrt() %>%
    matrix(size, t) %>%
    t()

  lim_i <- f1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- f1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  f1 <- as.data.frame(f1)
  lim_i <- as.data.frame(lim_i)
  lim_s <- as.data.frame(lim_s)

  names(f1) <- names_var
  names(lim_i) <- names_var
  names(lim_s) <- names_var

  max_value <- calcula_max(f1 - min(f1))[[3]] + min(f1)
  min_value <- -calcula_max(-(f1 - max(f1)))[[3]] + max(f1)

  if (max_value - min_value < 1e-2) {
    center <- (max_value + min_value) / 2
    max_value <- center + 0.01
    max_value <- center - 0.01
  }

  f1$time <- c(1:dim(f1)[1])
  lim_i$time <- c(1:dim(lim_i)[1])
  lim_s$time <- c(1:dim(lim_s)[1])

  f1 <- f1[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename("media" = "value")
  lim_i <- lim_i[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename("lim_i" = "value")
  lim_s <- lim_s[-c(1:cut_off), ] %>%
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

  plt <- ggplot(plot_data, aes_string(x = "time", fill = "fill_name", color = "color_name")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous("Time") +
    scale_color_manual("", values = color_list, na.value = NA) +
    scale_fill_manual("", values = fill_list, na.value = NA) +
    labs(title = paste0(pred, " (", ifelse(smooth, "smoothed", "only filtered"), ")")) +
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
