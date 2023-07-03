#### Default Method ####

#' Gamma outcome for kDGLM models
#'
#' Creates an outcome with gamma distribution with the chosen parameters (can only specify 2).
#'
#' @param phi character or numeric: Name of the linear predictor associated with the shape parameter of the gamma distribution. If numeric, this parameter is treated as known and equal to the value passed. If a character, the parameter is treated as unknown and equal to the exponential of the associated linear predictor. It cannot be specified with alpha.
#' @param mu character: Name of the linear predictor associated with the mean parameter of the gamma distribution. The parameter is treated as unknown and equal to the exponential of the associated linear predictor.
#' @param alpha character: Name of the linear predictor associated with the shape parameter of the gamma distribution. The parameter is treated as unknown and equal to the exponential of the associated linear predictor. It cannot be specified with phi.
#' @param beta character: Name of the linear predictor associated with the rate parameter of the gamma distribution. The parameter is treated as unknown and equal to the exponential of the associated linear predictor. It cannot be specified with sigma.
#' @param sigma character: Name of the linear predictor associated with the scale parameter of the gamma distribution. The parameter is treated as unknown and equal to the exponential of the associated linear predictor. It cannot be specified with beta.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @details
#'
#' For evaluating the posterior parameters, we use the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' @seealso \code{\link{fit_model}}
#' @family {auxiliary functions for a creating outcomes}
#'
#' @examples
#'
#' # Gamma case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' mu <- exp((sin(w * 1:T / T) + 2))
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
#' @references
#'    \insertAllCited{}
Gamma <- function(phi = NA, mu = NA, alpha = NA, beta = NA, sigma = NA, outcome, offset = outcome**0, alt_method = FALSE) {
  t <- length(outcome)
  r <- 1

  if (is.numeric(phi)) {
    k <- 1
    if (any(!is.na(c(alpha, beta, sigma)))) {
      stop("Error: When phi is known only mu can be estimated.")
    }
    pred_names <- c(mu)
    convert_mat_default <- convert_mat_canom <- diag(1)
    convert_canom_flag <- FALSE
    distr <- list(
      conj_prior = convert_Gamma_Normal,
      conj_post = convert_Normal_Gamma,
      update = update_Gamma,
      smoother = generic_smoother,
      calc_pred = gamma_pred,
      apply_offset = function(ft, Qt, offset) {
        t <- if.null(dim(ft)[2], 1)
        offset <- matrix(offset, t, r)

        list("ft" = ft + log(t(offset)), "Qt" = Qt)
      },
      link_function = log, inv_link_function = exp,
      param_names = c("alpha", "beta")
    )
    if (alt_method) {
      distr$update <- update_Gamma_alt
    }

    parms <- list(phi = phi)
  } else {
    warning("The estimation of the shape parameter phi is still under development. Results are not reliable.")
    k <- 2
    flags <- !is.na(c(phi, mu, alpha, beta, sigma))
    if (sum(flags) < 2) {
      stop("Error: Parameters not fully specified. You must specified exactly 2 non reduntant parameters.")
    }
    if (sum(flags) > 2) {
      stop("Error: Parameters over specified. You must specified exactly 2 non reduntant parameters.")
    }
    if (flags[4] & flags[5]) {
      stop("Error: Scale specified in more than one value.")
    }
    convert_mat_default <- matrix(c(1, 0, 0, 1, 1, 0, 1, -1, -2, 2), 2, 5)[, flags]
    convert_mat_canom <- solve(convert_mat_default)
    convert_canom_flag <- !all(flags[c(1, 2)])
    parms <- list()
    pred_names <- c(phi, mu, alpha, beta, sigma)[flags]
    names(pred_names) <- c("Shape (phi)", "Mean (mu)", "Shape (alpha)", "Rate (beta)", "Scale (sigma)")[flags]

    distr <- list(
      conj_prior = convert_FGamma_Normal,
      conj_post = convert_Normal_FGamma,
      update = update_FGamma,
      smoother = generic_smoother,
      calc_pred = Fgamma_pred,
      apply_offset = function(ft, Qt, offset) {
        list("ft" = ft + matrix(c(0, log(offset)), 2, dim(ft)[2]), "Qt" = Qt)
      },
      link_function = log, inv_link_function = exp,
      param_names = c("n", "k", "tau", "theta")
    )

    if (alt_method) {
      distr$conj_prior <- format_ft
      distr$conj_post <- format_param
      distr$update <- update_FGamma_alt
      distr$calc_pred <- Fgamma_pred_alt
      distr$param_names <- generic_param_names(k)
    }
  }

  distr$pred_names <- pred_names
  distr$r <- r
  distr$k <- k
  distr$l <- k
  distr$t <- t
  distr$offset <- matrix(offset, t, r)
  distr$outcome <- matrix(outcome, t, r)
  distr$convert_canom_flag <- convert_canom_flag
  distr$convert_mat_canom <- convert_mat_canom
  distr$convert_mat_default <- convert_mat_default
  distr$parms <- parms
  distr$name <- "Gamma"

  class(distr) <- "dlm_distr"
  distr$alt_method <- alt_method
  return(distr)
}

##### Gamma with unknown shape and mean #####

#' system_full_gamma
#'
#' Evaluate the compatibilizing equation for the full gamma model \insertCite{@see ArtigokParametrico;textual}{kDGLM}.
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @importFrom cubature cubintegrate
#'
#' @return A vector with the values of the system \insertCite{@see ArtigokParametrico;textual}{kDGLM}.
#' @keywords internal
system_full_gamma <- function(x, parms) {
  n <- exp(x) # exp(x[1])
  tau <- (parms$Hq1 - n) / parms$Hq2
  theta <- log(tau) - (1 + 5 / n) / (2 * parms$Hq1)

  a <- (n + 5) / 2
  b <- n(log(tau) - theta)

  # print((parms$Hq3 + parms$Hq4))
  if (n < 50) {
    # Densidade marginal aproximada de alpha (uso opcional).
    # f_densi=function(x){dgamma(a,b))}
    # c_val=1
    # Densidade marginal exata de phi.
    # f_densi_raw <- function(x) {
    #   exp(n * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - n * lgamma(x + 1) - (n * x + 1) * log(x * tau))
    # }
    f_densi_raw <- function(x) {
      exp(-n * lgamma(x) + x * n * theta + lgamma(n * x - 1) + log(x) - n * x * log(n * tau))
    }
    lim_sup <- Inf
    c_val <- cubintegrate(f_densi_raw, (1 + 1e-1) / n, lim_sup, nVec = 200)$integral
    # print("--------------------")
    # print(n)
    # print(c_val)
    f_densi <- function(x) {
      f_densi_raw(x) / c_val
    }
    # print('a')
    f <- function(x) {
      -x * (digamma(x * n - 1) - log(x * n * tau)) * f_densi(x)
    }
    Hp3 <- cubintegrate(f, 1 / n + 1e-2, lim_sup, nVec = 200)$integral

    # print('b')
    f <- function(x) {
      (-x * log(x) + lgamma(x)) * f_densi(x)
    }
    Hp4 <- cubintegrate(f, 1 / n + 1e-2, lim_sup, nVec = 200)$integral

    # print('c')
    # f <- function(x) {
    #   (x * digamma(x * n + 1)  - lgamma(x)- x*log(tau)) * f_densi(x)
    # }
    # Hp5 <- cubintegrate(f, 0, Inf, nVec = 200)$integral
    # print('sd')
    Hp5 <- Hp3 + Hp4


    # f <- function(x) {
    #   (lgamma(x)-x*digamma(n*x+1)) * f_densi(x)
    # }
    # Hp5 <- integrate(f, 0, Inf)$value+parms$Hq1*log(tau)
  } else {
    c_val <- 1
    # Hp3 <- log(tau / n) * a / b - 1 / n + b / (12 * (n**2) * (a - 1))
    # Hp4 <- a / b + 0.5 * (digamma(a) - log(b)) - b / (12 * (a - 1)) - 11 / 12
    # Hp5=Hp3+Hp4
    Hp5 <- a * log(tau) / b + 1 / n + 1 + b / (12 * (n**2) * (a - 1))
  }

  f_all <- c(
    (Hp5 - (parms$Hq3 + parms$Hq4))
  )
  # print(f_all)
  return(f_all)
}

a <- log(7e-18)

system_full_gamma2 <- function(x, parms) {
  n <- exp(x[1]) # exp(x[1])

  tau <- exp(x[2])
  theta <- exp(x[3])

  phi_proxy <- -3 + 3 * sqrt(1 - 4 * log(theta / tau) / 3)
  mu_proxy <- tau

  c <- n * (phi_proxy * log(phi_proxy / mu_proxy) -
    lgamma(phi_proxy) +
    log(theta) * phi_proxy - tau * phi_proxy / mu_proxy)


  f_densi_raw <- function(x) {
    phi <- x[1, ]
    mu <- x[2, ]

    val1 <- phi
    val2 <- phi / mu
    val3 <- phi * log(phi / mu) - lgamma(phi)

    l.fx <- n * (val3 + log(theta) * val1 - tau * val2) - c
    # a=max(l.fx,a)
    fx <- exp(l.fx - a) %>% as.matrix()

    rbind(
      fx,
      val1 * fx,
      val2 * fx,
      val3 * fx
    )
  }
  lim_sup <- Inf
  vals <- cubintegrate(f_densi_raw, c(0, 0), c(lim_sup, lim_sup), fDim = 4, nVec = 200)$integral

  Hp1 <- vals[2] / vals[1]
  Hp2 <- vals[3] / vals[1]
  Hp3 <- vals[4] / vals[1]
  # if (all(!is.nan(x))) {
  #   print(x)
  #   print(vals)
  # }
  return(c(
    parms$Hq1 - Hp1,
    parms$Hq2 - Hp2,
    parms$Hq3 + parms$Hq4 - Hp3
  ))
}

#' convert_FGamma_Normal
#'
#' Calculate the parameters of the conjugated prior that best approximates the given log-Normal distribution.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the conjugated prior.
#'
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this function.
#'
#' @importFrom rootSolve multiroot
#' @importFrom stats dlnorm
#'
#' @return The parameters of the conjugated distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
convert_FGamma_Normal <- function(ft, Qt, parms) {
  # s <- exp(ft[2, ] + 1)
  s <- 1
  f1 <- ft[1, ]
  f2 <- ft[2, ] - log(s)
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  Hq1 <- exp(f1 + q1 / 2)
  Hq2 <- exp(f1 - f2 + (q1 + q2 - 2 * q12) / 2)
  Hq3 <- (f2 + q12) * Hq1

  Hq4 <- cubintegrate(function(x) {
    (-x * log(x) + lgamma(x)) * dlnorm(x, f1, sqrt(q1))
  }, 0, Inf, nVec = 200)$integral

  parms <- list(
    "Hq1" = Hq1,
    "Hq2" = Hq2,
    "Hq3" = Hq3,
    "Hq4" = Hq4
  )

  # ss1 <- multiroot(f = system_full_gamma, start = c(0), parms = parms, maxiter = 2000)
  # ss1 <- multiroot(function(x){trigamma((exp(x)+5)/2)-q1},0)
  ss1 <- multiroot(f = system_full_gamma2, start = c(
    0,
    f2,
    -3 * ((exp(f1) / 3 + 1)**2) / 4 + 3 / 4 + log(tau)
  ), parms = parms, maxiter = 2000)


  x <- as.numeric(ss1$root)
  # n <- exp(x) # exp(x[1])
  # n=max(2/q1-5,1/Hq1+1e-2)

  # Calculando tau e theta dado n e k
  # tau <- (n * Hq1 - 1) / Hq2
  # theta <- n * log(tau / n) - (n + 5) / (2 * Hq1)
  # tau <- tau * s
  # theta <- theta + n * log(s)
  n <- exp(x[1])
  tau <- exp(x[2])
  theta <- exp(x[3])
  return(list("n" = n, "k" = n, "tau" = tau, "theta" = theta))
}

#' convert_Normal_FGamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given conjugated distribution.
#' The approximation is the best in the sense that it minimizes the KL divergence from the conjugated distribution to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the conjugated distribution (n, tau, theta).
#' @param parms list: A list of extra known parameters of the distribution. Not used in this function.
#'
#' @importFrom cubature cubintegrate
#'
#' @return The parameters of the Normal distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
convert_Normal_FGamma <- function(conj_prior, parms) {
  n <- conj_prior$n
  tau <- conj_prior$tau
  theta <- conj_prior$theta

  f_densi_raw <- function(x) {
    phi <- x[1, ]
    mu <- x[2, ]

    val1 <- phi
    val2 <- phi / mu
    val3 <- phi * log(phi / mu) - lgamma(phi)

    l.fx <- n * (val3 + log(theta) * val1 - tau * val2)
    # a=max(l.fx,a)
    fx <- exp(l.fx - a) %>% as.matrix()

    rbind(
      fx,
      log(phi) * fx,
      log(mu) * fx,
      (log(phi)**2) * fx,
      (log(mu)**2) * fx,
      log(phi) * log(mu) * fx,
    )
  }
  vals <- cubintegrate(f_densi_raw, c(0, 0), c(Inf, Inf), nVec = 200, fdim = 6)$integral

  ft <- matrix(c(vals[2], vals[3]), 2, 1) / vals[1]
  Qt <- matrix(c(vals[4], vals[6], vals[6], vals[5]), 2, 2) / vals[1] - ft %*% t(ft)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_FGamma
#'
#' Calculate posterior parameter for the conjugated distribution, assuming that the observed values came from a Gamma model with the shape (phi) and mean (mu) parameters having a conjugated prior.
#'
#' @param conj_prior list: A vector containing the parameters of the conjugated distribution (n, tau, theta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containing the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
update_FGamma <- function(conj_prior, ft, Qt, y, parms) {
  n0 <- conj_prior$n
  k0 <- conj_prior$k
  tau0 <- conj_prior$tau
  theta0 <- conj_prior$theta

  a <- (k0 + 1) / 2
  b <- (n0 - k0 + n0 * log(tau0 / n0) - theta0)

  n1 <- n0 + 1
  k1 <- k0 + 1
  tau1 <- tau0 + y
  theta1 <- theta0 + log(y)

  a <- (k1 + 1) / 2
  b <- (n1 - k1 + n1 * log(tau1 / n1) - theta1)

  return(list("n" = n1, "k" = k1, "tau" = tau1, "theta" = theta1))
}

#' Fgamma_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the distribution of the linear predictor.
#' The data is assumed to have Gamma distribution with unknown shape phi and unknown mean having log-Normal distribution.
#' In this scenario, the marginal distribution of the data is obtained via Monte Carlo.
#'
#' @param conj_param List or data.frame: The parameters of the distribution of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
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
#' @importFrom stats rgamma var quantile
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
Fgamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }
  r <- 1
  t <- length(conj_param$k)
  k <- 2

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
    }
    N <- 5000

    for (i in 1:t) {
      n <- conj_param$n[i]
      k <- conj_param$k[i]
      tau <- conj_param$tau[i]
      theta <- conj_param$theta[i]


      a <- (k + 1) / 2
      b <- (n - k + n * log(tau / n) - theta)

      alpha_i <- rgamma(N, a, b)
      mu_i <- 1 / rgamma(N, n * alpha_i + 1, alpha_i * tau)

      sample_y <- rgamma(N, alpha_i, alpha_i / mu_i)
      if (pred.flag) {
        pred[, i] <- mean(sample_y)
        var.pred[, , i] <- var(sample_y)
        icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
        icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
      }
      if (like.flag) {
        l.phi <- log(alpha_i)
        phi <- alpha_i
        l.mu <- log(mu_i)
        mu <- mu_i

        log.like.list <- phi * (l.phi - l.mu) - lgamma(phi) + phi * log(outcome[, i]) - phi * outcome[, i] / mu
        max.log.like <- max(log.like.list)
        like.list <- exp(log.like.list - max.log.like)
        log.like[i] <- log(mean(like.list)) + max.log.like
      }
    }
  }

  outcome_list <- list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
  return(outcome_list)
}

##### Gamma with known shape but unknown mean #####

#' convert_Gamma_Normal
#'
#' Calculate the parameters of the Inverse-Gamma that best approximates the given log-Normal distribution.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Inverse-Gamma
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this function.
#'
#' @return The parameters of the conjugated distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
convert_Gamma_Normal <- function(ft, Qt, parms) {
  alpha <- 1 / (-3 + 3 * sqrt(1 + 2 * Qt / 3))
  beta <- alpha * exp(ft - Qt / 2)
  return(list("alpha" = alpha, "beta" = beta))
}

#' convert_Normal_Gamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribution.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribution. Not used in this function.
#'
#' @return The parameters of the Normal distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
convert_Normal_Gamma <- function(conj_prior, parms) {
  ft <- -digamma(conj_prior$alpha) + log(conj_prior$beta)
  Qt <- trigamma(conj_prior$alpha)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Gamma
#'
#' Calculate posterior parameter for the Inverse-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribution Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
update_Gamma <- function(conj_prior, ft, Qt, y, parms) {
  alpha <- conj_prior$alpha + parms$phi
  beta <- conj_prior$beta + y * parms$phi
  return(list("alpha" = alpha, "beta" = beta))
}

#' gamma_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the conjugated distribution of the linear predictor.
#' The data is assumed to have Gamma distribution with known shape parameter phi and it's mean having distribution Inverse-Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribution of the data is Beta prime with parameters phi, alpha, beta / phi.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distributions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
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
#' @importFrom extraDistr qbetapr dbetapr
#' @export
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
gamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  phi <- parms$phi

  alpha <- conj_param$alpha %>% t()
  beta <- conj_param$beta %>% t()
  pred <- NULL
  var.pred <- NULL
  icl.pred <- NULL
  icu.pred <- NULL
  log.like <- NULL

  if (pred.flag) {
    pred <- ifelse(alpha > 1,
      beta / (alpha - 1),
      NA
    )
    var.pred <- ifelse(alpha > 1,
      ((beta / phi)**2) * (phi * (phi + alpha - 1) / ((alpha - 2) * (alpha - 1)**2)),
      Inf
    )


    icl.pred <- qbetapr((1 - pred_cred) / 2, phi, alpha, beta / phi)
    icu.pred <- qbetapr(1 - (1 - pred_cred) / 2, phi, alpha, beta / phi)
  }
  if (like.flag) {
    log.like <- dbetapr(outcome, phi, alpha, beta / phi, log = TRUE)
  }
  return(list(
    "pred" = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  ))
}

#### Alternative Method ####

##### Gamma with unknown shape and mean #####

#' update_FGamma_alt
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Gamma model from which the shape and mean parameters have prior distribution in the log-Normal family.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
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
update_FGamma_alt <- function(conj_prior, ft, Qt, y, parms) {
  # ft=c(0,0)
  # Qt=diag(2)
  # Qt[2,2]=2
  # y=6.063814

  if (all(diag(Qt) <= 0.1)) {
    # print("hey!")
    f0 <- ft
    S0 <- ginv(Qt)

    d1.log.like <- function(x) {
      phi <- exp(x[1])
      mu <- exp(x[2])
      # print(x)

      c(
        (log(phi) + 1 - log(mu) - digamma(phi) + log(y) - y / mu) * phi,
        (-phi / mu + phi * y / (mu**2)) * mu
      ) - S0 %*% (x - f0)
    }

    d2.inv <- function(x) {
      phi <- exp(x[1])
      mu <- exp(x[2])

      cross <- (-1 + y / mu) * phi

      mat <- matrix(
        c(
          (log(phi) + 1 - log(mu) - digamma(phi) + log(y) - y / mu + 1 - phi * trigamma(phi)) * phi,
          cross,
          cross,
          (-1 + y / mu + 1 - 2 * y / mu) * phi
        ),
        2, 2
      ) - S0
      return(mat)
    }

    mean <- c(f0[1], log(y))

    tau <- -d2.inv(mean) - S0

    f_start <- ginv(tau + S0) %*% (tau %*% mean + S0 %*% f0)
    # f_start=c(0.2286328,1.1518655)

    mode <- f_root(d1.log.like, d2.inv, f_start)$root
    # print(d2.inv(f0))
    # mode <- rootSolve::multiroot(d1.log.like, f_start)$root
    H <- d2.inv(mode)
    S <- ginv(-H)
  } else {
    f <- function(x) {
      l.phi <- x[1, ]
      l.mu <- x[2, ]
      phi <- exp(x[1, ])
      mu <- exp(x[2, ])

      # l.phi=log(x[1,])
      # l.mu=log(x[2,])
      # phi=x[1,]
      # mu=x[2,]
      # print(y)
      # print(dgamma(y, phi, phi / mu, log = TRUE))

      # phi*(l.phi-l.mu)-lgamma(phi)+phi*log(y)-phi*y/mu

      # prob <- exp(dgamma(y, phi, phi / mu, log = TRUE) + dmvnorm(t(x), ft, Qt, logged = TRUE))
      prob <- exp(phi * (l.phi - l.mu) - lgamma(phi) + phi * log(y) - phi * y / mu + dmvnorm(t(x), ft, Qt, logged = TRUE))


      rbind(
        prob,
        l.phi * prob,
        l.mu * prob,
        (l.phi**2) * prob,
        (l.phi * l.mu) * prob,
        (l.mu * l.phi) * prob,
        (l.mu**2) * prob
      )
    }

    # val <- cubintegrate(f, c(exp(mu2-6*sqrt(sigma2))), c(exp(mu2+6*sqrt(sigma2))), fDim = 7, nVec = 1000)$integral
    val <- cubintegrate(f, c(-Inf, -Inf), c(Inf, Inf), fDim = 7, nVec = 1000)$integral
    mode <- matrix(val[2:3] / val[1], 2, 1)
    S <- matrix(val[4:7], 2, 2) / val[1] - mode %*% t(mode)
  }

  return(list("ft" = matrix(mode, 2, 1), "Qt" = S))
}

#' Fgamma_pred_alt
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the distribution of the linear predictor.
#' The data is assumed to have Gamma distribution with unknown shape and mean parameters having log-Normal distribution.
#' In this scenario, the marginal distribution of the data is obtained via Monte Carlo.
#'
#' @param conj_param List or data.frame: The parameters of the distribution of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#' }
#'
#'
#' @importFrom stats rnorm rgamma var quantile
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
Fgamma_pred_alt <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)

  norm_param <- format_param(conj_param)
  ft <- norm_param$ft
  Qt <- norm_param$Qt

  k <- 2
  t <- dim(ft)[2]
  r <- 1

  pred <- NULL
  var.pred <- NULL
  icl.pred <- NULL
  icu.pred <- NULL
  log.like <- NULL

  Qt <- array(Qt, c(k, k, t))

  if (pred.flag | like.flag) {
    if (pred.flag) {
      pred <- matrix(NA, r, t)
      var.pred <- array(NA, c(r, r, t))
      icl.pred <- matrix(NA, r, t)
      icu.pred <- matrix(NA, r, t)
    }
    if (like.flag) {
      outcome <- matrix(outcome, r, t)
      log.like <- rep(NA, t)
    }

    N <- 5000
    sample <- matrix(rnorm(k * N), N, k)
    for (i in 1:t) {
      ft_i <- sample %*% var_decomp(Qt[, , i]) + matrix(ft[, i], N, k, byrow = TRUE)
      sample_y <- rgamma(N, exp(ft_i[, 1]), exp(ft_i[, 1] - ft_i[, 2]))
      if (pred.flag) {
        pred[, i] <- mean(sample_y)
        var.pred[, , i] <- var(sample_y)
        icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
        icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
      }
      if (like.flag) {
        l.phi <- ft_i[1, ]
        phi <- exp(l.phi)
        l.mu <- ft_i[2, ]
        mu <- exp(l.mu)

        log.like.list <- phi * (l.phi - l.mu) - lgamma(phi) + phi * log(outcome[, i]) - phi * outcome[, i] / mu
        max.log.like <- max(log.like.list)
        like.list <- exp(log.like.list - max.log.like)
        log.like[i] <- log(mean(like.list)) + max.log.like
      }
    }
  }

  outcome_list <- list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
  return(outcome_list)
}

##### Gamma with known shape but unknown mean #####

#' update_Gamma_alt
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Gamma model from which the shape parameter is known and mean parameter have prior distribution in the log-Normal family.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @importFrom cubature cubintegrate
#' @importFrom stats dlnorm
#'
#' @details
#'
#' For evaluating the posterior parameters, we use a modified version of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the detail about the modification of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, see \insertCite{ArtigoAltMethod;textual}{kDGLM}.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
#'
#' @references
#'    \insertAllCited{}
update_Gamma_alt <- function(conj_prior, ft, Qt, y, parms) {
  f <- function(x) {
    phi <- parms$phi
    l.phi <- log(phi)
    l.mu <- ft_i[2, ]
    mu <- exp(l.mu)

    log.like.list <- phi * (l.phi - l.mu) - lgamma(phi) + phi * log(y) - phi * y / mu

    # prob <- exp(dgamma(y, parms$phi, parms$phi / x, log = TRUE) + dlnorm(x, ft, sqrt(Qt), log = TRUE))
    prob <- exp(phi * (-l.mu) - phi * y / mu + dlnorm(x, ft, sqrt(Qt), log = TRUE))

    rbind(
      prob,
      log(x) * prob,
      (log(x)**2) * prob
    )
  }

  val <- cubintegrate(f, c(0), c(Inf), fDim = 3, nVec = 1000)$integral
  ft <- matrix(val[2] / val[1], 1, 1)
  Qt <- matrix(val[3], 1, 1) / val[1] - ft %*% t(ft)

  return(list("ft" = ft, "Qt" = Qt))
}
