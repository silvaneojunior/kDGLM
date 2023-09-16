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
#' @param data vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as data.
#'
#' @return A object of the class dlm_distr
#' @importFrom stats dgamma
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
#' structure <- polynomial_block(mu = 1, D = 0.95)
#'
#' outcome <- Gamma(phi = 0.5, mu = "mu", data = cornWheat$corn.log.return[1:500]**2)
#' fitted.data <- fit_model(structure, outcome)
#' summary(fitted.data)
#' plot(fitted.data, plot.pkg = "base")
#'
#' @references
#'    \insertAllCited{}
Gamma <- function(phi = NA, mu = NA, alpha = NA, beta = NA, sigma = NA, data, offset = as.matrix(data)**0) {
  alt.method <- FALSE
  data <- as.matrix(data)

  # phi=deparse(substitute(phi))[[1]] |> check.expr()
  # mu=deparse(substitute(mu))[[1]] |> check.expr()
  # alpha=deparse(substitute(alpha))[[1]] |> check.expr()
  # beta=deparse(substitute(beta))[[1]] |> check.expr()
  # sigma=deparse(substitute(sigma))[[1]] |> check.expr()

  t <- length(data)
  r <- 1

  if (is.numeric(phi)) {
    k <- 1
    if (any(!is.na(c(alpha, beta, sigma)))) {
      stop("Error: When phi is known only mu can be estimated.")
    }
    pred.names <- c("Mean (mu)"=mu)
    convert.mat.default <- convert.mat.canom <- diag(1)
    convert.canom.flag <- FALSE
    distr <- list(
      conj_distr = convert_Gamma_Normal,
      norm_distr = convert_Normal_Gamma,
      update = update_Gamma,
      smoother = generic_smoother,
      calc_pred = gamma_pred,
      apply_offset = function(ft, Qt, offset) {
        t <- if.null(dim(ft)[2], 1)
        offset <- matrix(offset, t, r)

        list("ft" = ft + log(t(offset)), "Qt" = Qt)
      },
      link_function = log, inv_link_function = exp,
      log_like_function = function(x, param) {
        dgamma(x, phi, phi / param)
      },
      param.names = c("alpha", "beta")
    )
    if (alt.method) {
      distr$update <- update_Gamma_alt
    }

    parms <- list(phi = phi)
  } else {
    if (!alt.method) {
      stop("Error: The estimation of the shape parameter phi is still under development. See the nightly build in GitHub to use this functionality.")
    }
    warning("The estimation of the shape parameter phi is still under development. Results are not reliable.")
    k <- 2
    flags <- !is.na(c(phi, mu, alpha, beta, sigma))
    if (sum(flags) < 2) {
      stop("Error: Parameters not fully specified. You must specified exactly 2 non reduntant parameters.")
    }
    if (sum(flags) > 2) {
      stop("Error: Parameters over specified. You must specified exactly 2 non reduntant parameters.")
    }
    if (all(flags[4:5])) {
      stop("Error: Scale specified in more than one variable.")
    }
    if (all(flags[c(1,3)])) {
      stop("Error: Shape specified in more than one variable.")
    }
    convert.mat.default <- matrix(c(1, 0, 0, 1, 1, 0, 1, -1, -2, 2), 2, 5)[, flags]
    convert.mat.canom <- solve(convert.mat.default)
    convert.canom.flag <- !all(flags[c(1, 2)])
    parms <- list()
    pred.names <- c(phi, mu, alpha, beta, sigma)[flags]
    names(pred.names) <- c("Shape (phi)", "Mean (mu)", "Shape (alpha)", "Rate (beta)", "Scale (sigma)")[flags]

    distr <- list(
      # conj_distr = convert_FGamma_Normal,
      # norm_distr = convert_Normal_FGamma,
      # update = update_FGamma,
      smoother = generic_smoother,
      # calc_pred = Fgamma_pred,
      apply_offset = function(ft, Qt, offset) {
        list("ft" = ft + matrix(c(0, log(offset)), 2, dim(ft)[2]), "Qt" = Qt)
      },
      link_function = log, inv_link_function = exp,
      log_like_function = function(x, param) {
        dgamma(x, param[1], param[1] / param[2], log = TRUE)
      },
      param.names = c("n", "k", "tau", "theta")
    )

    if (alt.method) {
      distr$conj_distr <- format_ft
      distr$norm_distr <- format_param
      distr$update <- update_FGamma_alt
      distr$calc_pred <- Fgamma_pred_alt
      distr$param.names <- generic_param_names(k)
    }
  }

  distr$pred.names <- pred.names
  distr$r <- r
  distr$k <- k
  distr$l <- k
  distr$t <- t
  distr$offset <- matrix(offset, t, r)
  distr$data <- matrix(data, t, r)
  distr$convert.canom.flag <- convert.canom.flag
  distr$convert.mat.canom <- convert.mat.canom
  distr$convert.mat.default <- convert.mat.default
  distr$parms <- parms
  distr$name <- "Gamma"

  class(distr) <- "dlm_distr"
  distr$alt.method <- alt.method
  return(distr)
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
#' @param conj.param list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribution. Not used in this function.
#'
#' @return The parameters of the Normal distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
convert_Normal_Gamma <- function(conj.param, parms) {
  ft <- -digamma(conj.param$alpha) + log(conj.param$beta)
  Qt <- trigamma(conj.param$alpha)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Gamma
#'
#' Calculate posterior parameter for the Inverse-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribution Inverse-Gamma.
#'
#' @param conj.param list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
update_Gamma <- function(conj.param, ft, Qt, y, parms) {
  alpha <- conj.param$alpha + parms$phi
  beta <- conj.param$beta + y * parms$phi
  return(list("alpha" = alpha, "beta" = beta))
}

#' gamma_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the conjugated distribution of the linear predictor.
#' The data is assumed to have Gamma distribution with known shape parameter phi and its mean having distribution Inverse-Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribution of the data is Beta prime with parameters phi, alpha, beta / phi.
#'
#' @param conj.param List or data.frame: The parameters of the conjugated distributions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' @param pred.cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred.cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred.cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item log.like vector: the The log likelihood for the outcome given the conjugated parameters.
#' }
#'
#' @importFrom extraDistr qbetapr dbetapr
#' @export
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with known shape}
gamma_pred <- function(conj.param, outcome = NULL, parms = list(), pred.cred = 0.95) {
  pred.flag <- !is.na(pred.cred)
  like.flag <- !is.null(outcome)
  if (!like.flag && !pred.flag) {
    return(list())
  }

  phi <- parms$phi

  alpha <- conj.param$alpha |> t()
  beta <- conj.param$beta |> t()
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

    icl.pred <- qbetapr((1 - pred.cred) / 2, phi, alpha, beta / phi)
    icu.pred <- qbetapr(1 - (1 - pred.cred) / 2, phi, alpha, beta / phi)
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
#' @param conj.param list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @importFrom cubature cubintegrate
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
update_FGamma_alt <- function(conj.param, ft, Qt, y, parms) {
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

    f.start <- ginv(tau + S0) %*% (tau %*% mean + S0 %*% f0)
    # f.start=c(0.2286328,1.1518655)

    mode <- f_root(d1.log.like, d2.inv, f.start)$root
    # print(d2.inv(f0))
    # mode <- rootSolve::multiroot(d1.log.like, f.start)$root
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
      prob <- exp(phi * (l.phi - l.mu) - lgamma(phi) + phi * log(y) - phi * y / mu + dmvnorm(t(x), ft, Qt))


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
#' @param conj.param List or data.frame: The parameters of the distribution of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' @param pred.cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred.cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred.cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#' }
#'
#'
#' @importFrom stats rgamma var quantile
#' @importFrom Rfast matrnorm
#' @keywords internal
#' @family {auxiliary functions for a Gamma outcome with unknowned shape}
Fgamma_pred_alt <- function(conj.param, outcome = NULL, parms = list(), pred.cred = 0.95) {
  pred.flag <- any(!is.na(pred.cred))
  like.flag <- any(!is.null(outcome))

  norm.param <- format_param(conj.param)
  ft <- norm.param$ft
  Qt <- norm.param$Qt

  k <- 2
  t <- dim(ft)[2]
  r <- 1

  pred <- NULL
  var.pred <- NULL
  icl.pred <- NULL
  icu.pred <- NULL
  log.like <- NULL

  Qt <- array(Qt, c(k, k, t))

  if (pred.flag || like.flag) {
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
    sample <- matrnorm(N, k)
    for (i in seq_len(t)) {
      ft.i <- rmvnorm(N, ft[, i], Qt[, , i])
      sample.y <- rgamma(N, exp(ft.i[, 1]), exp(ft.i[, 1] - ft.i[, 2]))
      if (pred.flag) {
        pred[, i] <- mean(sample.y)
        var.pred[, , i] <- var(sample.y)
        icl.pred[, i] <- quantile(sample.y, (1 - pred.cred) / 2)
        icu.pred[, i] <- quantile(sample.y, 1 - (1 - pred.cred) / 2)
      }
      if (like.flag) {
        l.phi <- ft.i[1, ]
        phi <- exp(l.phi)
        l.mu <- ft.i[2, ]
        mu <- exp(l.mu)

        log.like.list <- phi * (l.phi - l.mu) - lgamma(phi) + phi * log(outcome[, i]) - phi * outcome[, i] / mu
        max.log.like <- max(log.like.list)
        like.list <- exp(log.like.list - max.log.like)
        log.like[i] <- log(mean(like.list)) + max.log.like
      }
    }
  }

  outcome.list <- list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
  return(outcome.list)
}

##### Gamma with known shape but unknown mean #####

#' update_Gamma_alt
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Gamma model from which the shape parameter is known and mean parameter have prior distribution in the log-Normal family.
#'
#' @param conj.param list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
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
update_Gamma_alt <- function(conj.param, ft, Qt, y, parms) {
  f <- function(x) {
    phi <- parms$phi
    l.phi <- log(phi)
    l.mu <- ft
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
