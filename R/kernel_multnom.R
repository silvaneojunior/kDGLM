#### Default Method ####

#' Multinom outcome for kDGLM models
#'
#' Creates an outcome with Multinomial distribution with the chosen parameters.
#'
#' @param p character: a vector with the name of the linear predictor associated with the probability of each category (except the base one, which is assumed to be the last).
#' @param data vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as data.
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
#' structure <- (
#'   polynomial_block(p = 1, order = 2, D = 0.95) +
#'     harmonic_block(p = 1, period = 12, D = 0.975) +
#'     noise_block(p=1,R1=0.1)+
#'     regression_block(p = as.Date('2013-09-1')) # Vaccine was introduced in September of 2013
#' ) * 4
#'
#' outcome <- Multinom(p = structure$pred.names, data = chickenPox[, c(2, 3, 4, 6, 5)])
#' fitted.data <- fit_model(structure, outcome)
#' summary(fitted.data)
#' plot(fitted.data, plot.pkg = "base")
#'
#' @seealso \code{\link{fit_model}}
#' @family {auxiliary functions for a creating outcomes}
#'
#' @references
#'    \insertAllCited{}
Multinom <- function(p, data, offset = as.matrix(data)**0) {
  alt.method <- FALSE
  data <- as.matrix(data)

  # p=deparse(substitute(p))[[1]] |> check.expr()

  t <- dim(data)[1]
  r <- dim(data)[2]
  k <- dim(data)[2] - 1
  if (length(p) != k) {
    stop(paste0("Error: Incorrect number of parameter, expected ", k, " got ", length(p), "."))
  }
  convert.mat.default <- convert.mat.canom <- diag(k)
  parms <- list()
  # pred.names=p
  names(p) <- paste0("Odds for category ", 1:k, " (p.", 1:k, ")")
  distr <- list(
    pred.names = p,
    r = r,
    k = k,
    l = r,
    t = t,
    offset = matrix(offset, t, r),
    data = matrix(data, t, r),
    convert.mat.canom = convert.mat.canom,
    convert.mat.default = convert.mat.default,
    convert.canom.flag = FALSE,
    parms = parms,
    name = "Multinomial",
    conj_distr = convert_Multinom_Normal,
    norm_distr = convert_Normal_Multinom,
    update = update_Multinom,
    smoother = generic_smoother,
    calc_pred = multnom_pred,
    apply_offset = function(ft, Qt, offset) {
      t <- dim(ft)[2]
      offset <- matrix(offset, r, t)
      offset.class <- offset[-r, ]
      offset.ref <- matrix(offset[r, ], r - 1, t, byrow = TRUE)
      return(list("ft" = ft + log(offset.class / offset.ref), "Qt" = Qt))
    },
    link_function = function(x) {
      t <- length(x)
      return(log(x[-t, ] / x[t, ]))
    },
    inv_link_function = function(x) {
      y <- exp(x)
      x_last <- 1 / (1 + colSums(y))
      return(rbind(y * x_last, x_last))
    },
    log_like_function = function(x, param) {
      dmultinom(x, size = sum(x), prob = param, log = TRUE)
    },
    param.names = paste0("alpha_", 1:r)
  )
  class(distr) <- "dlm_distr"
  distr$alt.method <- alt.method

  if (alt.method) {
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
  calc.helper <- 1 + sum(exp(ft))
  r <- length(ft) + 1

  H <- exp(ft) %*% t(exp(ft)) / (calc.helper**2)
  diag(H) <- diag(H) - exp(ft) / calc.helper

  media.log <-
    -log(calc.helper) + sum(diag(0.5 * (H %*% Qt)))

  system_multinom <- function(x) {
    x <- exp(x)
    digamma_last <- digamma(x[r] - sum(x[-r]))
    digamma_vec <- digamma(x)

    f.all <- ft - digamma_vec[-r] + digamma_last
    last.guy <- media.log - digamma_last + digamma_vec[r]

    f.all <- c(f.all, last.guy)

    return(f.all)
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
#' @param conj.param list: A vector containing the concentration parameters of the Dirichlet.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the Normal distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
convert_Normal_Multinom <- function(conj.param, parms = list()) {
  alpha <- conj.param
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
#' @param conj.param list: A vector containing the concentration parameters of the Dirichlet.
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
update_Multinom <- function(conj.param, ft, Qt, y, parms = list()) {
  alpha <- conj.param + y
  return(alpha)
}

#' multnom_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the conjugated distribution of the linear predictor.
#' The data is assumed to have Multinomial distribution with known number of trial N and the probability vector having distribution Dirichlet with parameters alpha_i.
#' In this scenario, the marginal distribution of the data is Dirichlet-Multinomial with parameters N and alpha_i.
#'
#' @param conj.param List or data.frame: The parameters of the conjugated distributions of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time. The value passed is used to compute N.
#' @param parms List (optional): A list of extra parameters for the model. Not used in this function.
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
#' @importFrom Rfast data.frame.to_matrix Lgamma colCumSums
#' @keywords internal
#' @family {auxiliary functions for a Multinomial outcome}
multnom_pred <- function(conj.param, outcome, parms = list(), pred.cred = 0.95) {
  pred.flag <- !is.na(pred.cred)
  like.flag <- !is.null(outcome)
  if (!like.flag && !pred.flag) {
    return(list())
  }

  if (is.null(dim(conj.param))) {
    conj.param <- conj.param |>
      data.frame.to_matrix() |>
      t()
  }

  t <- nrow(conj.param)
  k <- ncol(conj.param) - 1
  r <- ncol(conj.param)

  pred <- NULL
  var.pred <- NULL
  icl.pred <- NULL
  icu.pred <- NULL
  log.like <- NULL
  if (pred.flag || like.flag) {
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

    for (t_i in seq_len(t)) {
      outcome_t <- outcome[t_i, ]
      N <- sum(outcome_t)
      N <- max(N, 1)

      alpha <- conj.param[t_i, ] |> as.numeric()
      alpha0 <- sum(alpha)

      const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

      if (pred.flag) {
        # alpha=c(5,4,3,9)*5
        # alpha0 <- sum(alpha)
        # N <- 2000
        # r=4
        # pred.cred=0.95
        # t_i=1
        # icl.pred=matrix(NA,r,1)
        # icu.pred=matrix(NA,r,1)
        # icl.pred2=matrix(NA,r,1)
        # icu.pred2=matrix(NA,r,1)
        # pred=matrix(NA,r,1)
        # var.pred=array(NA,c(r,r,1))
        # const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

        p <- alpha / alpha0

        pred[, t_i] <- N * p
        var.pred[, , t_i] <- diag(N * p * (1 - p) * (N + alpha0) / (alpha0 + 1))
        for (i in 2:r) {
          for (j in 1:(i - 1)) {
            var.pred[i, j, t_i] <- var.pred[j, i, t_i] <- -N * p[i] * p[j] * (N + alpha0) / (alpha0 + 1)
          }
        }

        x_mat <- matrix(0:N, N + 1, r)
        alpha_mat <- matrix(alpha, N + 1, r, byrow = TRUE)
        x_alpha_mat <- x_mat + alpha_mat

        lgamma_x_alpha_mat <- lgamma(x_alpha_mat)
        lgamma_neg_x_alpha_mat <- lgamma(N+alpha0-x_alpha_mat)

        lgamma_alpha_mat <- matrix(lgamma(alpha), N + 1, r, byrow = TRUE)
        lgamma_neg_alpha_mat <- matrix(lgamma(alpha0-alpha), N + 1, r, byrow = TRUE)

        lgamma_x_mat <- matrix(lgamma(seq_len(N+1)), N + 1, r)
        lgamma_neg_x_mat <- matrix(lgamma(N:0+1), N + 1, r)

        prob_mat <-
          lgamma_x_alpha_mat +
          lgamma_neg_x_alpha_mat +
          - lgamma_alpha_mat +
          - lgamma_neg_alpha_mat +
          - lgamma_x_mat +
          - lgamma_neg_x_mat

        prob_mat <- exp(const + prob_mat)
        probs_acum <- colCumSums(prob_mat)
        icl.pred[, t_i] <- colSums(probs_acum < ((1 - pred.cred) / 2))
        icu.pred[, t_i] <- colSums(probs_acum < (1 - (1 - pred.cred) / 2))
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
#' @param conj.param list: A vector containing the concentration parameters of the Dirichlet. Not used in the alternative method.
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
update_Multinom_alt <- function(conj.param, ft, Qt, y, parms = list()) {
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
    pre.mat <- diag(r - 1)
    diag(pre.mat) <- p[-r]

    mat <- n * (p[-r] %*% t(p[-r])) - n * pre.mat +
      -S0
    mat
  }

  # Calculating good initialization
  alpha0 <- sum(y + 0.01)
  mean <- digamma(y + 0.01) - digamma(alpha0)
  var <- diag(trigamma(y + 0.01)) - trigamma(alpha0)

  mini.A <- diag(length(ft))
  A <- cbind(mini.A, -1)
  A.t <- rbind(mini.A, -1)

  mean <- A %*% mean
  var <- A %*% var %*% A.t
  tau <- ginv(var)

  f.start <- ginv(tau + S0) %*% (tau %*% mean + S0 %*% f0)

  mode <- f_root(d1.log.like, d2.log.like, start = f.start)$root
  H <- d2.log.like(mode)
  S <- ginv(-H)
  return(list("ft" = matrix(mode, length(mode), 1), "Qt" = S))
}
