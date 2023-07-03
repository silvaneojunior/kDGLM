#' Summary for a fitted kDGLM model
#'
#' Prints a report for a fitted_dlm object.
#'
#' @param fitted_dlm A fitted_dlm object.
#' @param t Integer: The time index for the latent states.
#' @param lag Integer: The number of steps ahead used for the evaluating the latent variables. Use lag<0 for the smoothed distribution, If lag==0 for the filtered distribution and lag=h for the h-step-ahead prediction.
#' @param metric_cutoff Integer: The cutoff time index for the metric calculation. Values before that time will be ignored.
#' @param pred_cred numeric: The credibility interval to be used for the interval score.
#'
#' @export
#' @importFrom stats pnorm
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
#' report_dlm(fitted_data)
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
report_dlm <- function(fitted_dlm, t = fitted_dlm$t, lag = -1, metric_cutoff = floor(fitted_dlm$t / 10), pred_cred = 0.95) {
  r <- length(fitted_dlm$outcomes)
  k <- dim(fitted_dlm$mt)[1]
  T <- dim(fitted_dlm$mt)[2]
  predictions <- eval_past(fitted_dlm, T = metric_cutoff:T, lag = lag, pred_cred = pred_cred, eval_pred = TRUE)
  distr_like <- sum(predictions$log.like)
  distr_rae <- mean(predictions$rae)
  distr_mae <- mean(predictions$mae)
  distr_mse <- mean(predictions$mse)
  distr_interval_score <- sum(predictions$interval.score)

  distr_names <- lapply(fitted_dlm$outcomes, function(x) {
    x$name
  })
  distr_names_len <- max(sapply(names(distr_names), nchar))

  coef_label <- if (lag < 0) {
    "smoothed"
  } else if (lag == 0) {
    "filtered"
  } else if (lag == 1) {
    "one-step-ahead prediction"
  } else {
    paste0(lag, "-steps-ahead prediction")
  }
  coef_names <- rep(NA, k)
  for (name in names(fitted_dlm$var_names)) {
    name_len <- length(fitted_dlm$var_names[[name]])
    name_i <- name
    if (name_len > 1) {
      len_var <- length(fitted_dlm$var_names[[name]])
      len_char <- floor(log10(len_var)) + 1
      name_i <- paste0(name, "_", formatC(1:len_var, width = len_char, flag = "0"))
    }
    coef_names[fitted_dlm$var_names[[name]]] <- name_i
  }
  len_names <- max(sapply(as.character(coef_names), function(x) {
    nchar(x)
  }))
  coef_names <- format(coef_names, width = len_names, justify = "r")

  mean_coef <- predictions$mt[, t - metric_cutoff]
  var_mat <- predictions$Ct[, , t - metric_cutoff]
  if (length(var_mat) == 1) {
    std_coef <- sqrt(abs(var_mat))
  } else {
    std_coef <- sqrt(abs(diag(var_mat)))
  }
  t_coef <- mean_coef / std_coef
  p_val <- 2 * (1 - pnorm(abs(mean_coef) / std_coef))
  status <- rep(" ", length(coef_names))
  status[p_val <= 0.01] <- "."
  status[p_val <= 0.05] <- "*"
  status[p_val <= 0.01] <- "**"
  status[p_val <= 0.001] <- "***"
  mean_coef <- format(round(mean_coef, 5), width = 8, justify = "l")
  std_coef <- format(round(std_coef, 5), width = 10, justify = "l")
  t_coef <- format(round(t_coef, 5), width = 7, justify = "l")
  p_val_str <- ifelse(p_val < 0.001,
    format(p_val, digits = 3, justify = "l", scientific = TRUE),
    format(round(p_val, 3), width = 8, justify = "l", scientific = FALSE)
  )
  p_val_str <- ifelse(p_val < 1e-12,
    "  <1e-12",
    p_val_str
  )

  distr_like <- ifelse(abs(distr_like) < 0.00001,
    format(distr_like, digits = 4, width = 14, justify = "l", scientific = TRUE),
    format(round(distr_like, 5), width = 14, justify = "l", scientific = FALSE)
  )
  distr_rae <- ifelse(abs(distr_rae) < 0.00001,
    format(distr_rae, digits = 4, width = 21, justify = "l", scientific = TRUE),
    format(round(distr_rae, 5), width = 21, justify = "l", scientific = FALSE)
  )
  distr_mae <- ifelse(abs(distr_mae) < 0.00001,
    format(distr_mae, digits = 4, width = 17, justify = "l", scientific = TRUE),
    format(round(distr_mae, 5), width = 17, justify = "l", scientific = FALSE)
  )
  distr_mse <- ifelse(abs(distr_mse) < 0.00001,
    format(distr_mse, digits = 4, width = 20, justify = "l", scientific = TRUE),
    format(round(distr_mse, 5), width = 20, justify = "l", scientific = FALSE)
  )
  distr_interval_score <- ifelse(abs(distr_interval_score) < 0.00001,
    format(distr_interval_score, digits = 4, width = 16, justify = "l", scientific = TRUE),
    format(round(distr_interval_score, 5), width = 16, justify = "l", scientific = FALSE)
  )

  cat(paste0(
    "Fitted DGLM with ", length(fitted_dlm$outcomes), " outcomes.\n\n",
    "distributions:\n",
    paste0("    ", names(distr_names), ": ", distr_names, "\n", collapse = ""), "\n",
    "Coeficients (", coef_label, ") at time ", t, ":\n",
    paste(format(" ", width = len_names, justify = "l"), "Estimate", "Std. Error", " t value", "Pr(>|t|)"), "\n",
    paste(coef_names, mean_coef, std_coef, t_coef, p_val_str, status, "\n", collapse = ""),
    "---\n",
    "Signif. codes:  0 \xe2\x80\x98***\xe2\x80\x99 0.001 \xe2\x80\x98**\xe2\x80\x99 0.01 \xe2\x80\x98*\xe2\x80\x99 0.05 \xe2\x80\x98.\xe2\x80\x99 0.1 \xe2\x80\x98 \xe2\x80\x99 1\n\n",
    "---\n",
    format(" ", width = distr_names_len, justify = "l"), "  Pred. log-like  Mean abs. Error  Relative abs. Error  Mean Squared Error  Interval Score\n",
    paste0(format(names(distr_names), width = distr_names_len, justify = "l"), ": ", distr_like, distr_mae, distr_rae, distr_mse, distr_interval_score, "\n", collapse = ""),
    "---"
  ))
}

#' Summary for a kDGLM outcome
#'
#' Prints a report for a dlm_distr object.
#'
#' @param dlm_distr A dlm_distr object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for a creating outcomes}
report_distr <- function(dlm_distr) {
  cat(paste0(
    dlm_distr$name, " distribution.\n\nUnknown parameters: \n",
    paste0("    ", names(dlm_distr$pred_names), " - ", dlm_distr$pred_names, collapse = "\n"), "\n",
    if (length(dlm_distr$parms) > 0) {
      paste0("Known parameters: \n", paste0("    ", names(dlm_distr$parms), "=", dlm_distr$parms, collapse = "\n"), "\n")
    } else {
      ""
    },
    "\n",
    paste0("Serie length: ", dlm_distr$t, "\n"),
    paste0("Number of outcomes: ", dlm_distr$r, "\n"),
    paste0("Number of parameters: ", dlm_distr$k), "\n"
  ))
}

#' Summary for a kDGLM structure
#'
#' Prints a report for a dlm_block object.
#'
#' @param dlm_block A dlm_block object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for structural blocks}
report_block <- function(dlm_block) {
  block_names <- names(dlm_block$var_names)

  for (name in unique(block_names)) {
    count_name <- sum(block_names == name)
    if (count_name > 1) {
      len_char <- floor(log10(count_name)) + 1
      block_names[block_names == name] <- paste0(name, "_", formatC(1:count_name, width = len_char, flag = "0"))
    }
  }

  cat(paste0(
    dlm_block$type, " DLM block.",
    "\n",
    paste0("Latent variables: \n", paste0("    ", block_names, ": ", lapply(dlm_block$var_names, function(x) {
      paste0(names(x), collapse = ", ")
    }), " (", lapply(dlm_block$var_names, length), " variable(s))", collapse = "\n"), "\n"),
    "\n",
    paste0("Linear predictors: \n", paste0("    ", dlm_block$pred_names, collapse = "\n"), "\n"),
    "\n",
    paste0("Status: ", dlm_block$status, "\n"),
    paste0("Serie length: ", dlm_block$t, "\n"),
    paste0("Number of latent variables: ", dlm_block$n, "\n"),
    paste0("Number of linear predictors: ", dlm_block$k)
  ))
}

#' Summary for a searched_dlm object
#'
#' Prints a report for a searched_dlm object.
#'
#' @param searched_dlm A searched_dlm object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted_dlm class}
report_searched_dlm <- function(searched_dlm) {
  print(searched_dlm$search.data[1:5, ])
  report_dlm(searched_dlm$model)
}
