#' Visualizing a fitted kDGLM model
#'
#' Calculate the predictive mean and some quantile for the observed data and show a plot.
#'
#' @param model fitted_dlm: A fitted DGLM.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param lag Integer: A integer with the number of steps ahead should be used for prediction. If lag<0, the smoothed distribution is used and, if lag==0, the filtered distribuition is used.
#' @param plot_pkg String: A flag indicating if a plot should be produced. Should be one of 'auto', 'base', 'ggplot2' or 'plotly'.
#'
#' @return A list containing:
#' \itemize{
#'    \item data tibble object: A data frame containing the observations, predictions and credibility intervals at each time.
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#' }
#' @export
#' @importFrom grDevices rainbow
#' @import graphics
#' @import grDevices
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
#' show_fit(fitted_data, lag = -1)$plot
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
show_fit <- function(model, pred_cred = 0.95, lag = 1, plot_pkg = "auto") {
  if (plot_pkg == "auto") {
    plot_pkg <- if (requireNamespace("plotly", quietly = TRUE) & requireNamespace("ggplot2", quietly = TRUE)) {
      "plotly"
    } else if (requireNamespace("ggplot2", quietly = TRUE)) {
      "ggplot2"
    } else {
      "base"
    }
  }

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

  if (plot_pkg == "base" | !requireNamespace("ggplot2", quietly = TRUE)) {
    points <- paste0(colors, "55")
    fills <- paste0(colors, "33")
    names(colors) <- names(points) <- names(fills) <- series_names

    if (plot_pkg != "base") {
      warning("The ggplot2 package is required for ggplot2 and plotly plots and was not found. Falling back to R base plot functions.")
    }
    cur_height <- dev.size("cm")[2]
    count_spaces <- ceiling(n_colors / 4)
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
    plot(0, 0, type = "n", xlim = c(1, t_last), ylim = c(min_value, max_value), ylab = "$y_t$", xlab = "Time", main = title)
    for (serie in series_names) {
      plot_serie <- eval[eval$Serie == serie, ]
      points(plot_serie$Time[1:t_last], plot_serie$Observation[1:t_last],
        col = points[[serie]],
        pch = 16
      )
      lines(plot_serie$Time[1:t_last], plot_serie$Prediction[1:t_last], col = colors[[serie]])
      base_ribbon(plot_serie$Time, plot_serie$C.I.lower, plot_serie$C.I.upper,
        col = fills[[serie]], lty = 0
      )
    }

    par(mar = c(0, 0, 0, 0), cex = 1)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    legend(
      legend = c("Fit", "Obs."),
      col = c("black", "black"),
      lty = c(1, 0),
      seg.len = 0.6,
      pch = c(0, 16),
      pt.cex = c(0, 1),
      fill = c("#00000033", "#ffffff00"),
      border = "#ffffff00",
      cex = 0.75,
      x = "top", bty = "n", ncol = 2
    )
    par(mar = c(0, 0, 0, 0), cex = 1)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    legend(
      legend = series_names,
      col = colors,
      lty = rep(0, n_colors),
      pch = rep(22, n_colors),
      pt.cex = rep(2, n_colors),
      pt.bg = colors,
      x = 0.5, xjust = 0.5, y = 1, inset = 0, bty = "n",
      cex = 0.75,
      ncol = min(4, ceiling(n_colors / count_spaces))
    )
    par(mar = config$mar)
    plt <- NULL
  } else {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = eval, ggplot2::aes_string(x = "Time", fill = "Serie", ymin = "C.I.lower", ymax = "C.I.upper"), alpha = 0.25) +
      ggplot2::geom_line(data = eval, ggplot2::aes_string(x = "Time", color = "Serie", y = "Prediction", linetype = "'Fitted values'")) +
      ggplot2::geom_point(data = eval, ggplot2::aes_string(x = "Time", color = "Serie", y = "Observation", shape = '"Observation"'), alpha = 0.5) +
      ggplot2::scale_linetype_manual("", values = linetypes) +
      ggplot2::scale_shape_manual("", values = shapes) +
      ggplot2::scale_fill_manual("", na.value = NA, values = colors) +
      ggplot2::scale_color_manual("", na.value = NA, values = colors) +
      ggplot2::scale_y_continuous(name = "$y_t$") +
      ggplot2::scale_x_continuous("Time") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_bw() +
      ggplot2::coord_cartesian(ylim = c(min_value, max_value))
    for (name_i in names(model$outcomes)) {
      outcome_i <- model$outcomes[[name_i]]
      if (any(outcome_i$alt.flags == 1)) {
        plt <- plt +
          ggplot2::geom_vline(
            data = data.frame(xintercept = (1:t_last)[outcome_i$alt.flags == 1], linetype = "Detected changes", serie = name_i),
            ggplot2::aes_string(xintercept = "xintercept", linetype = "linetype")
          )
      }
    }
    if (plot_pkg == "plotly") {
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
  }
  return(list(data = eval, plot = plt))
}

#' Visualizing latent states in a fitted kDGLM model
#'
#' @param model fitted_dlm: A fitted DGLM model.
#' @param var Character: The name of the variables to plot (same value passed while creating the structure). Any variable whose name partially match this variable will be plotted.
#' @param lag Integer: A integer with the number of steps ahead should be used for evaluating the latent variables. If lag<0, the smoothed distribution is used and, if lag==0, the filtered distribution is used.
#' @param cutoff Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param plot_pkg String: A flag indicating if a plot should be produced. Should be one of 'auto', 'base', 'ggplot2' or 'plotly'.
#'
#' @return A list containing:
#' \itemize{
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#'    \item data data.frame: A data frame containing the data used in the plot.
#' }
#' @export
#'
#' @importFrom grDevices rainbow
#' @importFrom stats qnorm
#' @import graphics
#' @import grDevices
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
#' plot_lat_var(fitted_data)$plot
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
plot_lat_var <- function(model, var = "", lag = -1, cutoff = floor(model$t / 10), pred_cred = 0.95, plot_pkg = "auto") {
  if (plot_pkg == "auto") {
    plot_pkg <- if (requireNamespace("plotly", quietly = TRUE) & requireNamespace("ggplot2", quietly = TRUE)) {
      "plotly"
    } else if (requireNamespace("ggplot2", quietly = TRUE)) {
      "ggplot2"
    } else {
      "base"
    }
  }

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

  m1 <- param_distr$mt[indice, (cutoff + 1):t, drop = FALSE] |>
    t()
  std_mat <- param_distr$Ct[indice, indice, (cutoff + 1):t, drop = FALSE] |>
    apply(3, diag) |>
    sqrt() |>
    matrix(size, t - cutoff) |>
    t()

  lim_i <- m1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- m1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  max_value <- evaluate_max(m1 - min(m1))[[3]] + min(m1)
  min_value <- -evaluate_max(-(m1 - max(m1)))[[3]] + max(m1)

  if (max_value - min_value < 1e-2) {
    center <- (max_value + min_value) / 2
    max_value <- center + 0.01
    max_value <- center - 0.01
  }

  plot_data <- data.frame(
    Time = c((cutoff + 1):t),
    Label = as.factor(c(sapply(var_names, function(x) {
      rep(x, t - cutoff)
    }))),
    Mean = c(m1),
    C.I.lower = c(lim_i),
    C.I.upper = c(lim_s)
  )

  label <- paste0("C.I. (", pred_cred * 100 |> round(), "%)")

  var_names <- levels(plot_data$Label)
  n_var <- length(var_names)
  color_list <- rainbow(n_var, s = 0.6)
  names(color_list) <- var_names

  fill_list <- rainbow(n_var, s = 0.6)
  names(fill_list) <- var_names
  plot_data$fill_name <- plot_data$Label
  plot_data$color_name <- plot_data$Label

  title <- if (lag < 0) {
    "Smoothed estimation of latent variables"
  } else if (lag == 0) {
    "Filtered estimation of latent variables"
  } else if (lag == 1) {
    "One-step-ahead prediction for latent variables"
  } else {
    paste0(lag, "-steps-ahead prediction for latent variables")
  }

  if (plot_pkg == "base" | !requireNamespace("ggplot2", quietly = TRUE)) {
    points <- paste0(color_list, "55")
    fills <- paste0(color_list, "33")
    names(color_list) <- names(points) <- names(fills) <- var_names

    if (plot_pkg != "base") {
      warning("The ggplot2 package is required for ggplot2 and plotly plots and was not found. Falling back to R base plot functions.")
    }

    cur_height <- dev.size("cm")[2]
    count_spaces <- ceiling(n_var / 4)
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
    plot(0, 0, type = "n", xlim = c(cutoff, t), ylim = c(min_value, max_value), ylab = "$y_t$", xlab = "Time", main = title)
    for (var_name in var_names) {
      plot_serie <- plot_data[plot_data$Label == var_name, ]
      points(plot_serie$Time, plot_serie$Observation,
        col = points[[var_name]],
        pch = 16
      )
      lines(plot_serie$Time, plot_serie$Mean, col = color_list[[var_name]])
      base_ribbon(plot_serie$Time, plot_serie$C.I.lower, plot_serie$C.I.upper,
        col = fills[[var_name]], lty = 0
      )
    }
    lines(c(-1, t + 2), c(0, 0), lty = 2)

    par(mar = c(0, 0, 0, 0), cex = 1)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    legend(
      legend = c("Mean", label),
      col = c("black", "black"),
      lty = c(1, 0),
      seg.len = 0.6,
      pch = c(0, 22),
      pt.cex = c(0, 2),
      pt.bg = c("#00000000", "#00000033"),
      border = "#ffffff00",
      x = 0.5, xjust = 0.5, y = 1, bty = "n", cex = 0.75, horiz = TRUE
    )
    par(mar = c(0, 0, 0, 0), cex = 1)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    legend(
      legend = var_names,
      col = color_list,
      lty = rep(0, n_var),
      pch = rep(22, n_var),
      pt.cex = rep(2, n_var),
      pt.bg = color_list,
      x = 0.5, xjust = 0.5, y = 1, inset = 0, cex = 0.75, bty = "n",
      ncol = min(4, ceiling(n_var / count_spaces))
    )
    par(mar = config$mar)
    plt <- NULL
  } else {
    plt <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "Time", fill = "Label", color = "color_name")) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_x_continuous("Time") +
      ggplot2::scale_color_manual("", values = color_list, na.value = NA) +
      ggplot2::scale_fill_manual("", values = fill_list, na.value = NA) +
      ggplot2::labs(title = title) +
      ggplot2::scale_y_continuous("Parameter value") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "C.I.lower", ymax = "C.I.upper"), alpha = 0.25, color = NA) +
      ggplot2::geom_line(ggplot2::aes_string(y = "Mean")) +
      ggplot2::coord_cartesian(ylim = c(min_value, max_value))
    if (plot_pkg == "plotly") {
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
  }
  return(list(data = plot_data[, 1:5], plot = plt))
}

#' Visualizing linear predictors in a fitted kDGLM model
#'
#' @param model fitted_dlm: A fitted DGLM model.
#' @param pred Character: The name of the linear predictors to plot (same value passed while creating the structure). Any predictors whose name partially match this variable will be plotted.
#' @param lag Integer: A integer with the number of steps ahead should be used for evaluating the linear predictors. If lag<0, the smoothed distribution is used and, if lag==0, the filtered distribution is used.
#' @param cutoff Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not reliable.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param plot_pkg String: A flag indicating if a plot should be produced. Should be one of 'auto', 'base', 'ggplot2' or 'plotly'.
#'
#' @return A list containing:
#' \itemize{
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#'    \item data tibble: A data frame containing the data used in the plot.
#' }
#' @export
#' @importFrom grDevices rainbow
#' @importFrom stats qnorm
#' @import graphics
#' @import grDevices
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
#' plot_lin_pred(fitted_data)$plot
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
plot_lin_pred <- function(model, pred = "", lag = -1, cutoff = floor(model$t / 10), pred_cred = 0.95, plot_pkg = "auto") {
  if (plot_pkg == "auto") {
    plot_pkg <- if (requireNamespace("plotly", quietly = TRUE) & requireNamespace("ggplot2", quietly = TRUE)) {
      "plotly"
    } else if (requireNamespace("ggplot2", quietly = TRUE)) {
      "ggplot2"
    } else {
      "base"
    }
  }

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

  f1 <- param_distr$ft[indice, (cutoff + 1):t, drop = FALSE] |>
    t()
  std_mat <- param_distr$Qt[indice, indice, (cutoff + 1):t, drop = FALSE] |>
    apply(3, diag) |>
    sqrt() |>
    matrix(size, t - cutoff) |>
    t()

  lim_i <- f1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- f1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  max_value <- evaluate_max(f1 - min(f1))[[3]] + min(f1)
  min_value <- -evaluate_max(-(f1 - max(f1)))[[3]] + max(f1)

  if (max_value - min_value < 1e-2) {
    center <- (max_value + min_value) / 2
    max_value <- center + 0.01
    max_value <- center - 0.01
  }

  plot_data <- data.frame(
    Time = c((cutoff + 1):t),
    Label = as.factor(c(sapply(var_names, function(x) {
      rep(x, t - cutoff)
    }))),
    Mean = c(f1),
    C.I.lower = c(lim_i),
    C.I.upper = c(lim_s)
  )

  label <- paste0("C.I. (", pred_cred * 100 |> round(), "%)")

  n_var <- length(unique(plot_data$Label))
  color_list <- rainbow(n_var, s = 0.6)
  names(color_list) <- unique(plot_data$Label)

  fill_list <- rainbow(n_var, s = 0.6)
  names(fill_list) <- unique(plot_data$Label)
  plot_data$fill_name <- plot_data$Label
  plot_data$color_name <- plot_data$Label
  plot_data$IC_name <- plot_data$Label

  title <- if (lag < 0) {
    "Smoothed estimation of linear predictors"
  } else if (lag == 0) {
    "Filtered estimation of linear predictors"
  } else if (lag == 1) {
    "One-step-ahead prediction for linear predictors"
  } else {
    paste0(lag, "-steps-ahead prediction for linear predictors")
  }

  if (plot_pkg == "base" | !requireNamespace("ggplot2", quietly = TRUE)) {
    points <- paste0(color_list, "55")
    fills <- paste0(color_list, "33")
    names(color_list) <- names(points) <- names(fills) <- var_names

    if (plot_pkg != "base") {
      warning("The ggplot2 package is required for ggplot2 and plotly plots and was not found. Falling back to R base plot functions.")
    }

    cur_height <- dev.size("cm")[2]
    count_spaces <- ceiling(n_var / 4)
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
    plot(0, 0, type = "n", xlim = c(cutoff, t), ylim = c(min_value, max_value), ylab = "$y_t$", xlab = "Time", main = title)
    for (var_name in var_names) {
      plot_serie <- plot_data[plot_data$Label == var_name, ]
      points(plot_serie$Time, plot_serie$Observation,
        col = points[[var_name]],
        pch = 16
      )
      lines(plot_serie$Time, plot_serie$Mean, col = color_list[[var_name]])
      base_ribbon(plot_serie$Time, plot_serie$C.I.lower, plot_serie$C.I.upper,
        col = fills[[var_name]], lty = 0
      )
    }
    lines(c(-1, t + 2), c(0, 0), lty = 2)

    par(mar = c(0, 0, 0, 0), cex = 1)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    legend(
      legend = c("Mean", label),
      col = c("black", "black"),
      lty = c(1, 0),
      seg.len = 0.6,
      pch = c(0, 22),
      pt.cex = c(0, 2),
      pt.bg = c("#00000000", "#00000033"),
      border = "#ffffff00",
      x = 0.5, xjust = 0.5, y = 1, bty = "n", cex = 0.75, ncol = 2
    )
    par(mar = c(0, 0, 0, 0), cex = 1)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    legend(
      legend = var_names,
      col = color_list,
      lty = rep(0, n_var),
      pch = rep(22, n_var),
      pt.cex = rep(2, n_var),
      pt.bg = color_list,
      x = 0.5, xjust = 0.5, y = 1, inset = 0, cex = 0.75, bty = "n",
      ncol = min(4, ceiling(n_var / count_spaces))
    )
    par(mar = config$mar)
    plt <- NULL
  } else {
    plt <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "Time", fill = "Label", color = "color_name")) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_x_continuous("Time") +
      ggplot2::scale_color_manual("", values = color_list, na.value = NA) +
      ggplot2::scale_fill_manual("", values = fill_list, na.value = NA) +
      ggplot2::labs(title = title) +
      ggplot2::scale_y_continuous("predictor value") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "C.I.lower", ymax = "C.I.upper"), alpha = 0.25, color = NA) +
      ggplot2::geom_line(ggplot2::aes_string(y = "Mean")) +
      ggplot2::coord_cartesian(ylim = c(min_value, max_value))
    if (plot_pkg == "plotly") {
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
  }
  return(list(data = plot_data[, 1:5], plot = plt))
}
