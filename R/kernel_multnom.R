#### Default Method ####

#' Multinom outcome for kDGLM models
#'
#' Creates an outcome with Multinomial distribution with the chosen parameters.
#'
#' @param p character: a vector with the name of the linear predictor associated with the probality of each category (except the base one, which is assumed to be the last).
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#'
#' @importFrom stats dmultinom
#' @export
#'
#' @details
#'
#' For evaluating the posterior parameters, we use the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' @examples
#'
#' # Multinomial case
#' T <- 200
#' y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
#' y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
#' y3 <- rpois(T, exp(5))
#'
#' y <- cbind(y1, y2, y3)
#'
#' level1 <- polynomial_block(p1 = 1, order = 2)
#' level2 <- polynomial_block(p2 = 1, order = 2)
#' season <- harmonic_block(p2 = 1, period = 12)
#' outcome <- Multinom(p = c("p1", "p2"), outcome = y)
#'
#' fitted_data <- fit_model(level1, level2, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, lag = -1)$plot
#'
#' @seealso \code{\link{fit_model}}
#' @family {auxiliary functions for a creating outcomes}
#'
#' @references
#'    \insertAllCited{}
Multinom <- function(p, outcome, offset = outcome**0, alt_method = FALSE) {
  t <- dim(outcome)[1]
  r <- dim(outcome)[2]
  k <- dim(outcome)[2] - 1
  if (length(p) != k) {
    stop(paste0("Error: Incorrect number of parameter, expected ", k, " got ", length(p), "."))
  }
  convert_mat_default <- convert_mat_canom <- diag(k)
  parms <- list()
  # pred_names=p
  names(p) <- paste0("Odds for category ", c(1:k), " (p_", c(1:k), ")")
  distr <- list(
    pred_names = p,
    r = r,
    k = k,
    l = r,
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    convert_canom_flag = FALSE,
    parms = parms,
    name = "Multinomial",
    conj_prior = convert_Multinom_Normal,
    conj_post = convert_Normal_Multinom,
    update = update_Multinom,
    smoother = generic_smoother,
    calc_pred = multnom_pred,
    apply_offset = function(ft, Qt, offset) {
      t <- dim(ft)[2]
      offset <- matrix(offset, r, t)
      offset_class <- offset[-r, ]
      offset_ref <- matrix(offset[r, ], r - 1, t, byrow = TRUE)
      return(list("ft" = ft + log(offset_class / offset_ref), "Qt" = Qt))
    },
    link_function = function(x) {
      T <- length(x)
      return(log(x[-T, ] / x[T, ]))
    },
    inv_link_function = function(x) {
      y <- exp(x)
      x_last <- 1 / (1 + colSums(y))
      return(rbind(y * x_last, x_last))
    },
    param_names = paste0("alpha_", 1:r)
  )
  class(distr) <- "dlm_distr"
  distr$alt_method <- alt_method

  if (alt_method) {
    distr$update <- update_Multinom_alt
  }

  return(distr)
}

#' convert_Multinom_Normal
#'
#' Calculate the parameters of the Dirichlet that best approximates the given log-Normal distribution.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Dirichlet.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @importFrom Rfast transpose
#'
#' @return The parameters of the conjugated distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
convert_Multinom_Normal <- function(ft, Qt, parms = list()) {
  calc_helper <- 1 + sum(exp(ft))
  k <- length(ft)
  r <- k + 1

  H <- exp(ft) %*% t(exp(ft)) / (calc_helper**2)
  diag(H) <- diag(H) - exp(ft) / calc_helper

  media.log <-
    -log(calc_helper) + sum(diag(0.5 * (H %*% Qt)))

  system_multinom <- function(x) {
    x <- exp(x)
    digamma_last <- digamma(x[r] - sum(x[-r]))
    digamma_vec <- digamma(x)

    f_all <- ft - digamma_vec[-r] + digamma_last
    last_guy <- media.log - digamma_last + digamma_vec[r]

    f_all <- c(f_all, last_guy)

    return(f_all)
  }

  jacob_multinom <- function(x) {
    x <- exp(x)
    trigamma_last <- trigamma(x[r] - sum(x[-r]))
    trigamma_vec <- trigamma(x)

    jacob <- diag(trigamma_vec)
    jacob[-r, -r] <- -jacob[-r, -r] - trigamma_last
    jacob[r, -r] <- trigamma_last
    jacob[-r, r] <- trigamma_last
    jacob[r, r] <- jacob[r, r] - trigamma_last

    transpose(jacob * x)
  }

  p0 <- exp(c(ft, 0))
  p <- p0 / sum(p0)

  ss1 <- f_root(system_multinom,
    jacob_multinom,
    # start = log(c(rep(0.01, r-1), 0.01 * r)))
    start = log(c(p[-r], 1) * 1)
  )

  tau <- exp(as.numeric(ss1$root))

  alpha <- tau
  alpha[r] <- tau[r] - sum(tau[-r])
  return(alpha)
}

#' convert_Normal_Multinom
#'
#' Calculate the parameters of the log-Normal that best approximates the given Dirichlet distribution.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Dirichlet to the log-Normal
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the Normal distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
convert_Normal_Multinom <- function(conj_prior, parms = list()) {
  alpha <- conj_prior
  r <- length(alpha)
  k <- r - 1
  ft <- digamma(alpha[-r]) - digamma(alpha[r])
  Qt <- matrix(trigamma(alpha[r]), k, k)
  diag(Qt) <- trigamma(alpha[-r]) + trigamma(alpha[r])
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Multinom
#'
#' Calculate posterior parameter for the Dirichlet, assuming that the observed values came from a Multinomial model from which the number of trials is known and the prior distribution for the probabilities of each category have joint distribution Dirichlet.
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet.
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @family {auxiliary functions for a Multinomial outcome}
update_Multinom <- function(conj_prior, ft, Qt, y, parms = list()) {
  r <- length(y)
  alpha <- conj_prior + y
  return(alpha)
}

#' multnom_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the conjugated distribution of the linear predictor.
#' The data is assumed to have Multinomial distribution with known number of trial N and the probability vector having distribution Dirichlet with parameters alpha_i.
#' In this scenario, the marginal distribution of the data is Dirichlet-Multinomial with parameters N and alpha_i.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distributions of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time. The value passed is used to compute N.
#' @param parms List (optional): A list of extra parameters for the model. Not used in this function.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item log.like vector: the The log likelihood for the outcome given the conjugated parameters.
#' }
#'
#' @importFrom Rfast data.frame.to_matrix Lgamma colCumSums
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
multnom_pred <- function(conj_param, outcome, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  if (is.null(dim(conj_param))) {
    conj_param <- conj_param |>
      data.frame.to_matrix() |>
      t()
  }

  t <- nrow(conj_param)
  k <- ncol(conj_param) - 1
  r <- ncol(conj_param)

  pred <- NULL
  var.pred <- NULL
  icl.pred <- NULL
  icu.pred <- NULL
  log.like <- NULL
  if (pred.flag | like.flag) {
    if (pred.flag) {
      pred <- matrix(NA, r, t)
      var.pred <- array(NA, c(r, r, t))
      icl.pred <- matrix(NA, r, t)
      icu.pred <- matrix(NA, r, t)
    }
    if (like.flag) {
      outcome <- matrix(outcome, t, r)
      log.like <- rep(NA, t)
    } else {
      outcome <- matrix(1 / r, t, r)
    }

    for (t_i in 1:t) {
      outcome_t <- outcome[t_i, ]
      N <- sum(outcome_t)
      N <- max(N, 1)

      alpha <- conj_param[t_i, ] |> as.numeric()
      alpha0 <- sum(alpha)

      const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

      if (pred.flag) {
        # alpha=c(5,4,3,9)*5
        # alpha0 <- sum(alpha)
        # N <- 2000
        # r=4
        # pred_cred=0.95
        # t_i=1
        # icl.pred=matrix(NA,r,1)
        # icu.pred=matrix(NA,r,1)
        # icl.pred2=matrix(NA,r,1)
        # icu.pred2=matrix(NA,r,1)
        # pred=matrix(NA,r,1)
        # var.pred=array(NA,c(r,r,1))
        # const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

        p <- alpha / alpha0
        p_var <- p * (1 - p) / (alpha0 + 1)

        pred[, t_i] <- N * p
        var.pred[, , t_i] <- diag(N * p * (1 - p) * (N + alpha0) / (alpha0 + 1))
        for (i in 2:r) {
          for (j in 1:(i - 1)) {
            var.pred[i, j, t_i] <- var.pred[j, i, t_i] <- -N * p[i] * p[j] * (N + alpha0) / (alpha0 + 1)
          }
        }

        lgamma_vec <- if ((N + 1) > 50) {
          Lgamma
        } else {
          lgamma
        }

        lgamma_mat <- if ((N + 1) * r > 50) {
          Lgamma
        } else {
          lgamma
        }

        x_mat <- matrix(0:N, N + 1, r)
        alpha_mat <- matrix(alpha, N + 1, r, byrow = TRUE)

        lgamma_x_plus_mat <- matrix(lgamma_vec(0:N + 1), N + 1, r)
        lgamma_alpha_mat <- matrix(lgamma_vec(alpha), N + 1, r, byrow = TRUE)
        lgamma_N_x_plus_mat <- matrix(lgamma_vec(N:0 + 1), N + 1, r)
        lgamma_alpha0_alpha_mat <- matrix(lgamma_vec(alpha0 - alpha), N + 1, r, byrow = TRUE)

        x_alpha_mat <- x_mat + alpha_mat

        lgamma_x_alpha_mat <- x_alpha_mat
        lgamma_N_alpha0_alpha_mat <- N + alpha0 - x_alpha_mat

        lgamma_x_alpha_mat <- lgamma_mat(x_alpha_mat)
        lgamma_N_alpha0_alpha_mat <- lgamma_mat(N + alpha0 - x_alpha_mat)

        prob_mat <- matrix(NA, N + 1, r)

        prob_mat <-
          lgamma_x_alpha_mat -
          lgamma_x_plus_mat -
          lgamma_alpha_mat +
          lgamma_N_alpha0_alpha_mat -
          lgamma_N_x_plus_mat -
          lgamma_alpha0_alpha_mat
        prob_mat <- exp(const + prob_mat)
        probs_acum <- colCumSums(prob_mat)
        icl.pred[, t_i] <- colSums(probs_acum < ((1 - pred_cred) / 2))
        icu.pred[, t_i] <- colSums(probs_acum < (1 - (1 - pred_cred) / 2))
      }
      if (like.flag) {
        log.like[t_i] <- const + sum(lgamma(outcome_t + alpha) - lgamma(outcome_t + 1) - lgamma(alpha))
      }
    }
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

#### Alternative Method ####

#' update_Multinom_alt
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Multinomial model from which the logit probabilities have prior distribution in the log-Normal family.
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet. Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
#'
#' @details
#'
#' For evaluating the posterior parameters, we use a modified version of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For computational efficiency, we also use a Laplace approximations to obtain the first and second moments of the posterior \insertCite{@see @TierneyKadane1 and @TierneyKadane2 }{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the detail about the modification of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, see \insertCite{ArtigoAltMethod;textual}{kDGLM}.
#'
#' @references
#'    \insertAllCited{}
update_Multinom_alt <- function(conj_prior, ft, Qt, y, parms = list()) {
  f0 <- ft
  S0 <- ginv(Qt)
  n <- sum(y)
  r <- length(y)

  log.like <- function(x) {
    p0 <- c(exp(x), 1)
    p <- p0 / sum(p0)

    sum(y * log(p)) - 0.5 * t(x - f0) %*% S0 %*% (x - f0)
  }

  d1.log.like <- function(x) {
    p0 <- c(x, 0)
    p0 <- p0 - max(p0)
    p0 <- exp(p0)
    p <- p0 / sum(p0)

    y[-r] - n * p[-r] +
      -S0 %*% (x - f0)
  }

  d2.log.like <- function(x) {
    p0 <- c(x, 0)
    p0 <- p0 - max(p0)
    p0 <- exp(p0)
    p <- p0 / sum(p0)
    pre_mat <- diag(r - 1)
    diag(pre_mat) <- p[-r]

    mat <- n * (p[-r] %*% t(p[-r])) - n * pre_mat +
      -S0
    mat
  }

  # Calculating good initialization
  alpha0 <- sum(y + 0.01)
  mean <- digamma(y + 0.01) - digamma(alpha0)
  var <- diag(trigamma(y + 0.01)) - trigamma(alpha0)

  mini_A <- diag(length(ft))
  A <- cbind(mini_A, -1)
  A_t <- rbind(mini_A, -1)

  mean <- A %*% mean
  var <- A %*% var %*% A_t
  tau <- ginv(var)

  f_start <- ginv(tau + S0) %*% (tau %*% mean + S0 %*% f0)

  mode <- f_root(d1.log.like, d2.log.like, start = f_start)$root
  H <- d2.log.like(mode)
  S <- ginv(-H)
  return(list("ft" = matrix(mode, length(mode), 1), "Qt" = S))
}
