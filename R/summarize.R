#' Summary for a fitted kDGLM model
#'
#' Prints a report for a fitted.dlm object.
#'
#' @param fitted.dlm A fitted.dlm object.
#' @param t Integer: The time index for the latent states.
#' @param lag Integer: The number of steps ahead used for the evaluating the latent variables. Use lag<0 for the smoothed distribution, If lag==0 for the filtered distribution and lag=h for the h-step-ahead prediction.
#' @param metric.lag Integer: The number of steps ahead used for the evaluating the predictions used when calculating metrics. Use metric.lag<0 for the smoothed distribution, If metric.lag==0 for the filtered distribution and metric.lag=h for the h-step-ahead prediction.
#' @param metric.cutoff Integer: The cutoff time index for the metric calculation. Values before that time will be ignored.
#' @param pred.cred numeric: The credibility interval to be used for the interval score.
#'
#' @export
#' @importFrom stats pnorm
#'
#' @examples
#'
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
#' # summary(fitted.data)
#' # or
#' report_dlm(fitted.data)
#'
#' @family {auxiliary visualization functions for the fitted.dlm class}
report_dlm <- function(fitted.dlm, t = fitted.dlm$t, lag = -1, metric.lag = 1, metric.cutoff = floor(fitted.dlm$t / 10), pred.cred = 0.95) {
  o <- length(fitted.dlm$outcomes)
  k <- fitted.dlm$k
  T_len <- fitted.dlm$t
  predictions <- eval_past(fitted.dlm, eval_t = (metric.cutoff + 1):T_len, lag = metric.lag, pred.cred = pred.cred, eval.pred = TRUE)
  distr.like <- sum(predictions$log.like)
  distr.rae <- mean(predictions$rae)
  distr.mae <- mean(predictions$mae)
  distr.mse <- mean(predictions$mse)
  distr.mase <- mean(predictions$mase)
  distr.interval.score <- mean(predictions$interval.score)
  predictions <- eval_past(fitted.dlm, eval_t = seq_len(T_len), lag = min(lag, 0), pred.cred = pred.cred, eval.pred = FALSE)

  distr.names <- lapply(fitted.dlm$outcomes, function(x) {
    x$name
  })
  distr.names.len <- max(sapply(names(distr.names), nchar))

  coef.label <- if (lag < 0) {
    "smoothed"
  } else if (lag == 0) {
    "filtered"
  } else if (lag == 1) {
    "one-step-ahead prediction"
  } else {
    paste0(lag, "-steps-ahead prediction")
  }

  metric.label <- if (metric.lag < 0) {
    "Smoothed predictions"
  } else if (metric.lag == 0) {
    "Filtered predictions"
  } else if (metric.lag == 1) {
    "One-step-ahead prediction"
  } else {
    paste0(metric.lag, "-steps-ahead prediction")
  }
  coef.names <- rep(NA, k)
  for (name in names(fitted.dlm$var.names)) {
    name.len <- length(fitted.dlm$var.names[[name]])
    name_i <- name
    if (name.len > 1) {
      len.var <- length(fitted.dlm$var.names[[name]])
      len.char <- floor(log10(len.var)) + 1
      name_i <- paste0(name, ".", formatC(names(fitted.dlm$var.names[[name]]), width = len.char, flag = "0"))
    }
    coef.names[fitted.dlm$var.names[[name]]] <- name_i
  }
  len.names <- max(sapply(as.character(coef.names), function(x) {
    nchar(x)
  }))
  coef.names <- format(coef.names, width = len.names, justify = "r")

  mean.coef <- predictions$mt[, t]
  var.mat <- predictions$Ct[, , t]
  if (length(var.mat) == 1) {
    std.coef <- sqrt(abs(var.mat))
  } else {
    std.coef <- sqrt(abs(diag(var.mat)))
  }
  t.coef <- mean.coef / std.coef
  p.val <- 2 * (1 - pnorm(abs(mean.coef) / std.coef))
  status <- rep(" ", length(coef.names))
  status[p.val <= 0.01] <- "."
  status[p.val <= 0.05] <- "*"
  status[p.val <= 0.01] <- "**"
  status[p.val <= 0.001] <- "***"
  mean.coef <- format(round(mean.coef, 5), width = 8, justify = "l")
  std.coef <- format(round(std.coef, 5), width = 10, justify = "l")
  t.coef <- format(round(t.coef, 5), width = 7, justify = "l")
  p.val.str <- ifelse(p.val < 0.001,
    format(p.val, digits = 3, justify = "l", scientific = TRUE),
    format(round(p.val, 3), width = 8, justify = "l", scientific = FALSE)
  )
  p.val.str <- ifelse(p.val < 1e-12,
    "  <1e-12",
    p.val.str
  )

  distr.like <- ifelse(abs(distr.like) < 0.00001 | abs(distr.like) > 1e5,
    format(distr.like, digits = 4, width = 14, justify = "l", scientific = TRUE),
    format(round(distr.like, 5), width = 14, justify = "l", scientific = FALSE)
  )
  distr.rae <- ifelse(abs(distr.rae) < 0.00001 | abs(distr.rae) > 1e5,
    format(distr.rae, digits = 4, width = 21, justify = "l", scientific = TRUE),
    format(round(distr.rae, 5), width = 21, justify = "l", scientific = FALSE)
  )
  distr.mae <- ifelse(abs(distr.mae) < 0.00001 | abs(distr.mae) > 1e5,
    format(distr.mae, digits = 4, width = 17, justify = "l", scientific = TRUE),
    format(round(distr.mae, 5), width = 17, justify = "l", scientific = FALSE)
  )
  distr.mse <- ifelse(abs(distr.mse) < 0.00001 | abs(distr.mse) > 1e5,
    format(distr.mse, digits = 4, width = 20, justify = "l", scientific = TRUE),
    format(round(distr.mse, 5), width = 20, justify = "l", scientific = FALSE)
  )
  distr.mase <- ifelse(abs(distr.mase) < 0.00001 | abs(distr.mase) > 1e5,
    format(distr.mase, digits = 4, width = 24, justify = "l", scientific = TRUE),
    format(round(distr.mase, 5), width = 24, justify = "l", scientific = FALSE)
  )
  distr.interval.score <- ifelse(abs(distr.interval.score) < 0.00001 | abs(distr.interval.score) > 1e5,
    format(distr.interval.score, digits = 4, width = 16, justify = "l", scientific = TRUE),
    format(round(distr.interval.score, 5), width = 16, justify = "l", scientific = FALSE)
  )

  cat(paste0(
    "Fitted DGLM with ", length(fitted.dlm$outcomes), " outcomes.\n\n",
    "distributions:\n",
    paste0("    ", names(distr.names), ": ", distr.names, "\n", collapse = ""), "\n",
    "Coeficients (", coef.label, ") at time ", t, ":\n",
    paste(format(" ", width = len.names, justify = "l"), "Estimate", "Std. Error", "  t value", "Pr(>|t|)"), "\n",
    paste(coef.names, mean.coef, std.coef, t.coef, p.val.str, status, "\n", collapse = ""),
    "---\n",
    "Signif. codes:  0 \xe2\x80\x98***\xe2\x80\x99 0.001 \xe2\x80\x98**\xe2\x80\x99 0.01 \xe2\x80\x98*\xe2\x80\x99 0.05 \xe2\x80\x98.\xe2\x80\x99 0.1 \xe2\x80\x98 \xe2\x80\x99 1\n\n",
    "---\n",
    metric.label, "\n",
    format(" ", width = distr.names.len, justify = "l"), "  Pred. log-like  Mean Abs. Error  Mean Abs. Scaled Error  Relative abs. Error  Mean Squared Error  Interval Score\n",
    paste0(format(names(distr.names), width = distr.names.len, justify = "l"), ": ", distr.like, distr.mae, distr.mase, distr.rae, distr.mse, distr.interval.score, "\n", collapse = ""),
    "---"
  ))
}

#' Summary for a kDGLM outcome
#'
#' Prints a report for a dlm_distr object.
#'
#' @param dlm.distr A dlm_distr object.
#'
#' @export
#'
#' @keywords internal
#' @family {auxiliary functions for a creating outcomes}
#' @family {Reports for dlm_distr objects.}
report_distr <- function(dlm.distr) {
  cat(paste0(
    dlm.distr$name, " distribution.\n\nUnknown parameters: \n",
    paste0("    ", names(dlm.distr$pred.names), " - ", dlm.distr$pred.names, collapse = "\n"), "\n",
    if (length(dlm.distr$parms) > 0) {
      paste0("Known parameters: \n", paste0("    ", names(dlm.distr$parms), "=", dlm.distr$parms, collapse = "\n"), "\n")
    } else {
      ""
    },
    "\n",
    paste0("Serie length: ", dlm.distr$t, "\n"),
    paste0("Number of outcomes: ", dlm.distr$r, "\n"),
    paste0("Number of parameters: ", dlm.distr$k), "\n"
  ))
}

#' Summary for a kDGLM structure
#'
#' Prints a report for a dlm_block object.
#'
#' @param dlm.block A dlm_block object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary functions for structural blocks}
report_block <- function(dlm.block) {
  block.names <- names(dlm.block$var.names)

  for (name in unique(block.names)) {
    count.name <- sum(block.names == name)
    if (count.name > 1) {
      len.char <- floor(log10(count.name)) + 1
      block.names[block.names == name] <- paste0(name, ".", formatC(1:count.name, width = len.char, flag = "0"))
    }
  }

  cat(paste0(
    dlm.block$type, " DLM block.",
    "\n",
    paste0("Latent variables: \n", paste0("    ", block.names, ": ", lapply(dlm.block$var.names, function(x) {
      paste0(names(x), collapse = ", ")
    }), " (", lapply(dlm.block$var.names, length), " variable(s))", collapse = "\n"), "\n"),
    "\n",
    paste0("Linear predictors: \n", paste0("    ", dlm.block$pred.names, collapse = "\n"), "\n"),
    "\n",
    paste0("Status: ", dlm.block$status, "\n"),
    paste0("Serie length: ", dlm.block$t, "\n"),
    paste0(
      "Interventions at: ",
      paste0(
        lapply(
          dlm.block$interventions,
          function(x) {
            paste0(x$times, collapse = ", ")
          }
        ),
        collapse = ", "
      ), "\n"
    ),
    paste0("Number of latent variables: ", dlm.block$n, "\n"),
    paste0("Number of linear predictors: ", dlm.block$k)
  ))
}

#' Summary for a searched_dlm object
#'
#' Prints a report for a searched_dlm object.
#'
#' @param searched.dlm A searched_dlm object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted.dlm class}
report_searched_dlm <- function(searched.dlm) {
  print(searched.dlm$search.data[1:5, ])
  report_dlm(searched.dlm$model)
}
