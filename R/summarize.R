#' Search the best hyper parameters fo a kDGLM model
#'
#' @param ... dlm_block object: The structural blocks of the model. At least one block must be "undefined". See polynomial_block for more details.
#' @param outcomes List: The observed data. It should contain objects of the class dlm_distr.
#' @param search_grid List: A named list containing the possible values of each undefined hyper parameter.
#' @param condition Character: A character defining which combinations of undefined hyper parameter should be tested. See example for details.
#' @param smooth Bool: A flag indicating if the smoothed distribution for the latent variables should be calculated.
#' @param pred_cred Numeric: A number between 0 and 1 (not included) indicating the credibility interval for predictions. If not within the valid range of values, 0.95 will be used.
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
#' This is an auxiliary function to search for the best combination of hyper parameters among the possible variations specfied by the user.
#' This function simple evaluates all possible combinations and chooses the one that results in the best fitting.
#' To compare each model, we evaluate the predictive log likelyhood for the one-step ahead prediction. We also compute the Relative Absolute Error (RAE) and return it to the user. It is important to note that one should be careful with the RAE metric, as it only takes into account the ponctual prediction.
#'
#' @seealso \code{\link{fit_model}}
search_model <- function(..., outcomes, search_grid, condition = "TRUE", smooth = TRUE, pred_cred = 0.95) {
  if (is.null(pred_cred) | is.na(pred_cred)) {
    pred_cred <- 0.95
  } else {
    if (pred_cred <= 0 | pred_cred >= 1) {
      warning("Invalid value for credibility. Using 95% of credibility.")
      pred_cred <- 0.95
    }
  }
  structure <- block_merge(...)
  if (structure$status == "defined") {
    stop("Error: There is no hiper parameter to select. Did you forgot to label the hiper parameters?")
  }

  ref_strucuture <- structure
  var_length <- length(search_grid)
  search_data <- do.call(expand.grid, search_grid) %>%
    filter(eval(parse(text = condition)))
  search_data$log.like <- NA
  search_data$rae <- NA
  vals_names <- names(search_grid)

  dim_FF <- dim(structure$FF)
  dimnames_FF <- dimnames(structure$FF)

  dim_D <- dim(structure$D)
  dim_W <- dim(structure$W)
  dim_G <- dim(structure$G)

  dim_C0 <- dim(structure$C0)

  vals_nums <- dim(search_data)[1]
  vals_size <- dim(search_data)[2]
  best_structure <- NULL
  best_model <- NULL
  best_metric <- -Inf
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
      structure$FF[structure$FF == name] <- cur_param[[name]]
      structure$D[structure$D == name] <- cur_param[[name]]
      structure$W[structure$W == name] <- cur_param[[name]]
      structure$m0[structure$m0 == name] <- cur_param[[name]]
      structure$C0[structure$C0 == name] <- cur_param[[name]]
      structure$G[structure$G == name] <- cur_param[[name]]
    }
    if (any(is.na(as.numeric(if.na(structure$FF, 0)))) |
      any(is.na(as.numeric(if.na(structure$D, 0)))) |
      any(is.na(as.numeric(if.na(structure$W, 0)))) |
      any(is.na(as.numeric(if.na(structure$m0, 0)))) |
      any(is.na(as.numeric(if.na(structure$C0, 0)))) |
      any(is.na(as.numeric(if.na(structure$G, 0))))
    ) {
      stop("Error: not all unkown hiper parameter have values. Check the search grid to make sure every unkown hiper parameter has a range of values.")
    }
    structure$FF <- array(as.numeric(structure$FF), dim_FF, dimnames = dimnames_FF)
    structure$D <- array(as.numeric(structure$D), dim_D)
    structure$W <- array(as.numeric(structure$W), dim_W)
    structure$m0 <- as.numeric(structure$m0)
    structure$C0 <- matrix(as.numeric(structure$C0), dim_C0[1], dim_C0[2])
    structure$G <- array(as.numeric(structure$G), dim_G)

    if (any(if.na(structure$D, 0) < 0 | if.na(structure$D, 0) > 1)) {
      stop(paste0("Error: invalid value for D. Expected a real number between 0 and 1, got: ", paste(structure$D[if.na(structure$D, 0) < 0 | if.na(structure$D, 0) > 1], collapse = ", "), "."))
    }
    if (any(if.na(structure$W, 0) < 0)) {
      stop(paste0("Error: invalid value for W. Expected a non negative number, got: ", paste(structure$W[if.na(structure$W, 0) < 0], collapse = ", "), "."))
    }
    if (any(if.na(structure$C0, 0) < 0)) {
      stop(paste0("Error: invalid value for C0. Expected a non negative number, got: ", paste(structure$C0[if.na(structure$C0, 0) < 0], collapse = ", "), "."))
    }

    structure$status <- "defined"
    fitted_model <- fit_model(structure, outcomes = outcomes, pred_cred = 0.95, smooth = FALSE, p_monit = NA, c_monit = 1)
    log.like <- sum(as.numeric(lapply(fitted_model$outcomes, function(x) {
      sum(x$log.like)
    })), na.rm = TRUE)
    rae <- mean(as.numeric(lapply(fitted_model$outcomes, function(x) {
      out <- x$outcome
      pred <- t(x$pred)

      base <- ifelse(out == 0, 1, out)
      mean(abs((pred - out) / out))
    })), na.rm = TRUE)
    search_data[i, vals_size - 1] <- log.like
    search_data[i, vals_size] <- rae
    if (log.like > best_metric) {
      best_structure <- structure
      best_model <- fitted_model
      best_metric <- log.like
    }
  }
  if (smooth) {
    smoothed <- generic_smoother(best_model$mt, best_model$Ct, best_model$at, best_model$Rt, best_model$G, best_model$G_labs)
    best_model$mts <- smoothed$mts
    best_model$Cts <- smoothed$Cts
    best_model$smooth <- smooth
  }
  cat("\n")
  search_data <- search_data %>% arrange(-log.like, -rae)

  out.vals <- list(
    search.data = search_data,
    structure = best_structure,
    model = best_model
  )
  class(out.vals) <- "searched_dlm"

  return(out.vals)
}

#' Summary for a fitted kDGLM model
#'
#' Prints a report for a fitted_dlm object.
#'
#' @param fitted_dlm A fitted_dlm object.
#' @param t Integer: The time index for the latent states.
#' @param smooth Bool: A flag indicating if the smoothed distribution should be used for the latent states.
#' @param metric_cutoff Integer: The cutoff time index for the metric calculation. Values before that time will be ignored.
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
report_dlm <- function(fitted_dlm, t = fitted_dlm$t, smooth = fitted_dlm$smooth, metric_cutoff = round(fitted_dlm$t / 10)) {
  r <- length(fitted_dlm$outcomes)
  k <- dim(fitted_dlm$mt)[1]
  distr_names <- list()
  distr_like <- rep(NA, r)
  distr_rae <- rep(NA, r)
  for (outcome_index in 1:r) {
    outcome <- fitted_dlm$outcomes[[outcome_index]]
    distr_names[names(fitted_dlm$outcomes)[outcome_index]] <- outcome$name

    prediction <- outcome$calc_pred(outcome$conj_prior_param[-(1:metric_cutoff), ], outcome$outcome[-(1:metric_cutoff), ], parms = outcome$parms, pred_cred = 0.95)
    distr_like[outcome_index] <- sum(prediction$log.like, na.rm = TRUE)

    pred <- t(prediction$pred)
    out <- outcome$outcome[-(1:metric_cutoff), ]
    out_div <- ifelse(out == 0, 1, out)
    distr_rae[outcome_index] <- mean(abs((pred - out) / out_div), na.rm = TRUE)
  }
  distr_names_len <- max(sapply(names(distr_names), nchar))

  coef_label <- if (fitted_dlm$smooth & smooth) {
    "smoothed"
  } else {
    "filtered"
  }
  coef_mean_name <- if (fitted_dlm$smooth & smooth) {
    "mts"
  } else {
    "mt"
  }
  coef_var_name <- if (fitted_dlm$smooth & smooth) {
    "Cts"
  } else {
    "Ct"
  }
  coef_names <- rep(NA,k)
  for (name in names(fitted_dlm$names)) {
    name_len <- length(fitted_dlm$names[[name]])
    name_i=name
    if (name_len > 1) {
      name_i <- paste0(name, "_", 1:length(fitted_dlm$names[[name]]))
    }
    coef_names[fitted_dlm$names[[name]]]=name_i
  }
  len_names <- max(sapply(as.character(coef_names), function(x) {
    nchar(x)
  }))
  coef_names <- format(coef_names, width = len_names, justify = "r")

  mean_coef <- fitted_dlm[[coef_mean_name]][, t]
  var_mat <- fitted_dlm[[coef_var_name]][, , t]
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

  cat(paste0(
    "Fitted DGLM with ", length(fitted_dlm$outcomes), " outcomes.\n\n",
    "distributions:\n",
    paste0(names(distr_names), ": ", distr_names, "\n", collapse = ""), "\n",
    "Coeficients (", coef_label, ") at time ", t, ":\n",
    paste(format(" ", width = len_names, justify = "l"), "Estimate", "Std. Error", " t value", "Pr(>|t|)"), "\n",
    paste(coef_names, mean_coef, std_coef, t_coef, p_val_str, status, "\n", collapse = ""),
    "---\n",
    "Signif. codes:  0 \xe2\x80\x98***\xe2\x80\x99 0.001 \xe2\x80\x98**\xe2\x80\x99 0.01 \xe2\x80\x98*\xe2\x80\x99 0.05 \xe2\x80\x98.\xe2\x80\x99 0.1 \xe2\x80\x98 \xe2\x80\x99 1\n\n",
    "---\n",
    format(" ", width = distr_names_len, justify = "l"), "  Pred. log-like  Relative abs. Error\n",
    paste0(format(names(distr_names), width = distr_names_len, justify = "l"), ": ", distr_like, distr_rae, "\n", collapse = ""),
    "---"
  ))
}

#' Summary for a kDGLM outcome
#'
#' Prints a report for a dlm_distr object.
#'
#' @param dlm_distr A fitted_dlm object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted_dlm class}
report_distr <- function(dlm_distr) {
  cat(paste0(
    dlm_distr$name, " distribution.\n\nUnkown parameters:\n",
    paste0(names(dlm_distr$var_names), ": ", dlm_distr$var_names, "\n", collapse = ""),
    if (length(dlm_distr$parms) > 0) {
      paste0(names(dlm_distr$parms), ": ", dlm_distr$parms, "\n", collapse = "\n")
    } else {
      "\n"
    },
    paste0("Serie length: ", dlm_distr$t, "\n"),
    paste0("Number of outcomes: ", dlm_distr$r)
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
