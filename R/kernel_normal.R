#### Default Method ####

#' Normal outcome for kDGLM models
#'
#' Creates an outcome with Normal distribution with the chosen parameters (can only specify 2,).
#'
#' @param mu character: Name of the linear predictor associated with the mean parameter of the Normal distribution. The parameter is treated as unknowed and equal to the associated linear predictor.
#' @param Tau character or numeric: If Tau is a character, it is interpreted as the names of the linear preditors associated with the precisions parameter of the Normal distribution. If Tau is numeric, the precision is considered known and equal to the value of Tau, otherwise, the precison is considere unknowned and equal to the exponential of the linear predictor informed in Tau. If the outcome is a Multivariate Normal, then Tau must be a matrix and, if the precision is unknowned, the elements outside it's main diagonal are treated as the linear predictor associated with the correlation between the each coordinate of the outcome, otherwise Tau is treated as the precision matrix. The user cannot be specify Tau with Sigma or Sd.
#' @param Sigma character or numeric: If Sigma is a character, it is interpreted as the names of the linear preditors associated with the variance parameter of the Normal distribution. If Sigma is numeric, the variance is considered known and equal to the value of Sigma, otherwise, the variance is considere unknowned and equal to the exponential of the linear predictor informed in Sigma. If the outcome is a Multivariate Normal, then Sigma must be a matrix and, if the variance is unknowned, the elements outside it's main diagonal are treated as the linear predictor associated with the correlation between the each coordinate of the outcome, otherwise Sigma is treated as the covariance matrix. The user cannot be specify Sigma with Tau or Sd.
#' @param Sd character or numeric: If Sd is a character, it is interpreted as the names of the linear preditors associated with the standard deviation parameter of the Normal distribution. If Sd is numeric, the standard deviation is considered known and equal to the value of Sd, otherwise, the precison is considere unknowned and equal to the exponential of the linear predictor informed by in Sd. If the outcome is a Multivariate Normal, then Tau must be a matrix and the elements outside it's main diagonal are treated as the correlation (or the name of the linear predictor associated) between the each coordinate of the outcome. The user cannot be specify Sd with Tau or Sigma.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @importFrom Rfast dmvnorm upper_tri.assign lower_tri.assign is.symmetric
#' @export
#'
#' @examples
#'
#' # Normal case
#' T <- 200
#' mu <- rnorm(T, 0, 0.1)
#' data <- rnorm(T, cumsum(mu))
#'
#' level <- polynomial_block(
#'   mu = 1,
#'   D = 0.95
#' )
#' variance <- polynomial_block(
#'   sigma2 = 1
#' )
#'
#' # Known variance
#' outcome <- Normal(mu = "mu", Sigma = 1, outcome = data)
#'
#' fitted_data <- fit_model(level, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Unknown variance
#' outcome <- Normal(mu = "mu", Sigma = "sigma2", outcome = data)
#'
#' fitted_data <- fit_model(level, variance, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
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
#' @references
#'    \insertAllCited{}
Normal <- function(mu, Tau = NA, Sigma = NA, Sd = NA, outcome, offset = outcome**0, alt_method = FALSE) {
  t <- if.null(dim(outcome)[1], length(outcome))
  r <- if.null(dim(outcome)[2], 1)

  #### Consistancy check ####
  {
    if (length(mu) != r) {
      stop("Error: The number of means does not match the number of series.")
    }
    if ((any(!is.na(Tau)) + any(!is.na(Sigma)) + any(!is.na(Sd))) > 1) {
      stop("Error: Can only specify one of Tau, Sigma or Sd")
    }
    if ((any(!is.na(Tau)) + any(!is.na(Sigma)) + any(!is.na(Sd))) == 0) {
      stop("Error: Must specify one of Tau, Sigma or Sd")
    }
    Var <- as.matrix(if (any(!is.na(Tau))) {
      Tau
    } else if (any(!is.na(Sigma))) {
      Sigma
    } else {
      Sd
    })
    Var_name <- as.matrix(if (any(!is.na(Tau))) {
      "Precision"
    } else if (any(!is.na(Sigma))) {
      "Variance"
    } else {
      "Standard deviation"
    })

    lower.index <- matrix(1:(r**2), r, r)[lower.tri(Var)]
    upper.index <- t(matrix(1:(r**2), r, r))[lower.tri(Var)]
    lower.flag <- lower.tri(Var)
    upper.flag <- upper.tri(Var)

    ###### Dimensions check ######
    if (!all(is.na(Var))) {
      if (!all(dim(Var) == r)) {
        stop(paste0("Error: The ", Var_name, " size is not compatible with the number of series. Expected ", r, "x", r, ", got ", paste0(dim(Var), collapse = "x")))
      }
    }

    ###### Symmetry test #####
    flags_up <- is.na(Var[upper.flag]) | is.null(Var[upper.flag])
    Var[upper.flag] <- ifelse(flags_up, Var[lower.flag], Var[upper.flag])

    flags_down <- is.na(Var[lower.flag]) | is.null(Var[lower.flag])
    Var[lower.flag] <- ifelse(flags_down, Var[upper.flag], Var[lower.flag])

    if (any(Var[lower.flag & !is.na(Var)] != Var[upper.flag & !is.na(Var)])) {
      if (!is.symmetric(Var)) {
        stop(paste0("Error: ", Var_name, " matrix is not symmetric."))
      }
    }

    ###### Under specification test #####
    if (any(is.na(Var))) {
      stop(paste0("Error: Covariance is not fully specified."))
    }
  }
  #### Specifying atributes and methods ####
  if (all(sapply(Var, is.numeric))) {
    var_names <- c(mu)
    names(var_names) <- c("mu")
    k <- r

    if (Var_name == "Standard deviation") {
      warning("Covariance matrix specified by the standard deviation. Non diagonal values are intepreted as correlation.")
    }

    Var <- if (Var_name == "Precision") {
      ginv(Tau)
    } else if (Var_name == "Variance") {
      Sigma
    } else {
      Cor <- Sd
      Sd <- diag(r)
      diag(Sd) <- diag(Cor)
      diag(Cor) <- 1
      Var <- diag(Sd) %*% Cor %*% diag(Sd)
    }
    distr <- list(
      conj_prior = format_ft,
      conj_post = format_param,
      update = update_Normal,
      smoother = generic_smoother,
      calc_pred = normal_pred,
      apply_offset = function(ft, Qt, offset) {
        t <- dim(ft)[2]
        offset <- t(matrix(offset, t, r))
        Qt <- array(Qt, c(r, r, t))

        ft_off <- ft * offset
        Qt_off <- sapply(1:t, function(i) {
          diag(offset[, i]) %*% Qt[, , i] %*% diag(offset[, i])
        }, simplify = "array")
        if (t == 1) {
          Qt_off <- matrix(Qt_off, r, r)
        }
        return(list("ft" = ft_off, "Qt" = Qt_off))
      },
      link_function = function(x) {
        x
      },
      inv_link_function = function(x) {
        x
      },
      param_names = generic_param_names(r)
    )
    parms <- list(Sigma = Var)
    convert_mat_canom <- convert_mat_default <- diag(r)
    convert_canom_flag <- FALSE

    distr$alt_method <- FALSE
  } else {
    k <- r + r * (r + 1) / 2
    # warning("Covariance matrix is unknowned. BewareNon diagonal values are intepreted as correlation.")

    var_names <- c(mu, diag(Var), Var[lower.index])
    mu_index <- 1:r
    var_index <- mu_index + r
    cor_index <- (2 * r + 1):k

    if (r > 1) {
      names_mat <- matrix(paste0("rho_", matrix(1:r, r, r, byrow = TRUE), matrix(1:r, r, r)), r, r)
      diag(names_mat) <- paste0(if (Var_name == "Precision") {
        "tau"
      } else if (Var_name == "Variance") {
        "sigma2"
      } else {
        "sigma"
      }, "_", 1:r)
      names(var_names) <- c(paste0("mu_", 1:r), diag(names_mat), names_mat[lower.index])
    } else {
      names(var_names) <- c("mu", if (Var_name == "Precision") {
        "tau"
      } else if (Var_name == "Variance") {
        "sigma2"
      } else {
        "sigma"
      })
    }
    distr <- list(
      conj_prior = if (r > 1) {
        convert_multi_NG_Normal
      } else {
        convert_NG_Normal
      },
      conj_post = if (r > 1) {
        convert_multi_Normal_NG
      } else {
        convert_Normal_NG
      },
      update = if (r > 1) {
        update_multi_NG_correl
      } else {
        if (alt_method) {
          update_NG_alt
        } else {
          update_NG
        }
      },
      smoother = generic_smoother,
      calc_pred = multi_normal_gamma_pred,
      apply_offset = function(ft, Qt, offset) {
        ft[mu_index, ] <- ft[mu_index, ] * offset
        ft[var_index, ] <- ft[var_index, ] - log(offset)

        Qt[mu_index, ] <- Qt[mu_index, ] * offset
        Qt[, mu_index] <- Qt[, mu_index] * offset
        return(
          list("ft" = ft, "Qt" = Qt)
        )
      },
      link_function = if (r > 1) {
        function(x) {
          x[var_index, ] <- log(x[var_index, ])
          a <- x[cor_index, ]
          a <- ((a + 1) / 2)
          a <- log(a / (1 - a))
          x[cor_index, ] <- a
          return(x)
        }
        # function(x) {
        #   x[var_index, ] <- log(x[var_index, ])
        #   a=matrix(1,r,r)
        #   a=upper_tri.assign(a,x[cor_index, ])
        #   a=lower_tri.assign(a,x[cor_index, ])
        #   A=eigen(a)
        #   vals=log(A$values)
        #   Q=A$vectors
        #   a <- Q%*%diag(vals)%*%t(Q)
        #   # a_raw
        #   # for(i in 1:r){
        #   #   a[r,]= a[r,]/sqrt(a_raw[r,r])
        #   #   a[,r]= a[,r]/sqrt(a_raw[r,r])
        #   # }
        #   x[cor_index, ] <- lower_tri(a)
        #   return(x)
        # }
      } else {
        function(x) {
          x[var_index, ] <- log(x[var_index, ])
          return(x)
        }
      },
      inv_link_function = if (r > 1) {
        function(x) {
          x[var_index, ] <- exp(x[var_index, ])
          a <- x[cor_index, ]
          a <- 1 / (1 + exp(-a))
          a <- 2 * a - 1
          x[cor_index, ] <- a
          return(x)
        }
        # function(x) {
        #   x[var_index, ] <- exp(x[var_index, ])
        #   a=matrix(0,r,r)
        #   a=upper_tri.assign(a,x[cor_index, ])
        #   a=lower_tri.assign(a,x[cor_index, ])
        #   A=eigen(a)
        #   vals=exp(A$values)
        #   Q=A$vectors
        #   a <- Q%*%diag(vals)%*%t(Q)
        #   a_raw=a
        #   for(i in 1:r){
        #     a[r,]= a[r,]/sqrt(a_raw[r,r])
        #     a[,r]= a[,r]/sqrt(a_raw[r,r])
        #   }
        #   x[cor_index, ] <- lower_tri(a)
        #   return(x)
        # }
      } else {
        function(x) {
          x[var_index, ] <- exp(x[var_index, ])
          return(x)
        }
      },
      param_names = c("mu0", "c0", "alpha0", "beta0")
    )
    if (r > 1) {
      distr$param_names <- paste(rep(distr$param_names, r), sort(rep(1:r, 4)), sep = "_")
    }

    parms <- list(
      alt_method = alt_method,
      mu_index = mu_index,
      var_index = var_index,
      cor_index = cor_index,
      upper.index = upper.index,
      lower.index = lower.index
    )
    convert_mat_canom <- convert_mat_default <- diag(k)
    convert_canom_flag <- FALSE
    if (Var_name == "Standard deviation") {
      diag(convert_mat_canom)[var_index] <- -2
      diag(convert_mat_default)[var_index] <- -0.5
      convert_canom_flag <- TRUE
    } else if (Var_name == "Variance") {
      diag(convert_mat_canom)[var_index] <- -1
      diag(convert_mat_default)[var_index] <- -1
      convert_canom_flag <- TRUE
    }

    distr$alt_method <- alt_method | r > 1
  }

  distr$var_names <- var_names
  distr$r <- r
  distr$k <- k
  distr$l <- k
  distr$t <- t
  distr$offset <- matrix(offset, t, r)
  distr$outcome <- matrix(outcome, t, r)
  distr$convert_mat_canom <- convert_mat_canom
  distr$convert_mat_default <- convert_mat_default
  distr$convert_canom_flag <- convert_canom_flag
  distr$parms <- parms
  distr$name <- "Normal"

  class(distr) <- "dlm_distr"

  return(distr)
}

#### Default Method ####
##### Normal with known variance #####

#' update_Normal
#'
#' Calculate posterior parameter for the Normal, assuming that the observed values came from a Normal model from which the covariance is known and the prior distribution for the mean vector have Normal distribution
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Normal.
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the covariance matrix parameter (Sigma) for the observational Normal model.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Normal outcome}
update_Normal <- function(conj_prior, ft, Qt, y, parms) {
  Sigma <- parms$Sigma
  if (length(Sigma) == 1) {
    null.flag <- Sigma == 0
  } else {
    null.flag <- all(diag(Sigma) == 0)
  }

  if (null.flag) {
    Qt <- Sigma * 0
    ft <- y
  } else {
    Tau0 <- ginv(Qt)
    Tau1 <- ginv(parms$Sigma)

    Qt <- ginv(Tau0 + Tau1)
    ft <- Qt %*% (Tau0 %*% ft + Tau1 %*% y)
  }
  return(do.call(c, list(ft, Qt)))
}

#' normal_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the conjugated distribution of the linear predictor.
#' The data is assumed to have Normal distribution with known variance and it's mean having distribution Normal.
#' In this scenario, the marginal distribution of the data is also Normal.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distributions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the observational covariance matrix, Sigma.
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
#' @importFrom Rfast dmvnorm data.frame.to_matrix
#' @importFrom stats qnorm
#' @keywords internal
#' @family {auxiliary functions for a Normal outcome}
normal_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  if (is.null(dim(conj_param))) {
    conj_param <- conj_param %>% matrix(1, length(conj_param))
  }
  r <- dim(conj_param)[2]
  r <- (-1 + sqrt(1 + 4 * r)) / 2
  t <- dim(conj_param)[1]

  pred <- t(data.frame.to_matrix(conj_param[, 1:r,drop=FALSE]))
  var.pred <- conj_param[(r + 1):(r * r + r)] %>%
    t() %>%
    array(c(r, r, t))
  icl.pred <- matrix(NA, t, r)
  icu.pred <- matrix(NA, t, r)
  log.like <- rep(NA, t)
  if (like.flag) {
    outcome <- matrix(outcome, t, r)
  }
  for (i in 1:t) {
    mu <- pred[, i]
    sigma2 <- (var.pred[, , i] + parms$Sigma) %>% matrix(r, r)
    if (pred.flag) {
      var.pred[, , i] <- sigma2
      # if (length(sigma2) > 1) {
      #   sigma2 <- diag(sigma2)
      # }
      icl.pred[i, ] <- rep(qnorm((1 - pred_cred) / 2), r) * sqrt(diag(sigma2)) + mu
      icu.pred[i, ] <- rep(qnorm(1 - (1 - pred_cred) / 2), r) * sqrt(diag(sigma2)) + mu
    }
    if (like.flag) {
      log.like[i] <- dmvnorm(outcome[i, ], mu, sigma2, logged = TRUE)
    }
  }
  if (!pred.flag) {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }
  if (!like.flag) {
    log.like <- NULL
  }

  list(
    "pred" = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

##### Normal with unknown variance #####


convert_NG_Normal <- function(ft, Qt, parms = list()) {
  ft <- matrix(ft, 2, 1)
  mu0 <- ft[1, ] + Qt[1, 2]
  c0 <- exp(-ft[2, ] - Qt[2, 2] / 2) / (Qt[1, 1] + 1e-10)
  helper <- -3 + 3 * sqrt(1 + 2 * Qt[2, 2] / 3)
  # helper=Qt[2,2]
  alpha <- 1 / helper
  beta <- alpha * exp(-ft[2, ] - Qt[2, 2] / 2)
  return(list("mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta))
}

convert_Normal_NG <- function(conj_distr, parms = list()) {
  # q12 <- conj_distr$mu0*(digamma(conj_distr$alpha)-log(conj_distr$beta))/(1-log(conj_distr$alpha)+log(conj_distr$beta))
  f1 <- conj_distr$mu0 #-q12
  f2 <- digamma(conj_distr$alpha) - log(conj_distr$beta + 1e-40)
  # f2 <- log(conj_distr$alpha)-log(conj_distr$beta)
  q1 <- conj_distr$beta / (conj_distr$c0 * conj_distr$alpha)
  q2 <- trigamma(conj_distr$alpha)
  q12 <- 0

  ft <- c(f1, f2)
  Qt <- matrix(c(q1, q12, q12, q2), byrow = F, ncol = 2)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_NG
#'
#' Calculate posterior parameter for the Normal-Gamma, assuming that the observed values came from a Normal model from which the prior distribution for the mean and the precision have joint distribution Normal-Gamma
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
update_NG <- function(conj_prior, ft, Qt, y, parms = list()) {
  mu0 <- (conj_prior$c0 * conj_prior$mu0 + y) / (conj_prior$c0 + 1)
  c0 <- conj_prior$c0 + 1
  alpha <- conj_prior$alpha + 0.5
  beta <- conj_prior$beta + 0.5 * conj_prior$c0 * ((conj_prior$mu0 - y)**2) / (conj_prior$c0 + 1)
  return(list("mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta))
}

##### Alternative method #####

#' update_NG_alt
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Normal model from which the mean and the log precision have prior distribution in the Normal family.
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @importFrom stats dnorm dlnorm
#' @keywords internal
#' @family {auxiliary functions for a Normal outcome}
#'
#' @details
#'
#' For evaluating the posterior parameters, we use a modified version of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the detail about the modification of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, see \insertCite{ArtigoAltMethod;textual}{kDGLM}.
#'
#' @references
#'    \insertAllCited{}
update_NG_alt <- function(conj_prior, ft, Qt, y, parms = list()) {
  # ft = c(-0.291027, -2.451121)
  # Qt = matrix(c(1.728035e-01, -1.199923e-03, -1.199923e-03,  2.461883e-05),2,2)
  # y=-1.234201
  # S0 <- ginv(Qt)
  # print(diag(Qt))

  if (Qt[2, 2] > 1e-2) {
    # if(flag_method){
    mu1 <- ft[1]
    mu2 <- ft[2]
    sigma1 <- Qt[1, 1]
    sigma2 <- Qt[2, 2]
    rho <- Qt[1, 2] / sqrt(sigma1 * sigma2)

    f <- function(x2) {
      mu1_cond <- mu1 + rho * sqrt(sigma1 / sigma2) * (log(x2) - mu2)
      sigma1_cond <- (1 - rho**2) * sigma1

      prec_obs <- x2
      prec_mu <- 1 / sigma1_cond
      mu1_cond_update <- (prec_obs * y + prec_mu * mu1_cond) / (prec_obs + prec_mu)
      sigma1_cond_update <- 1 / (prec_obs + prec_mu)

      prob <- exp(dnorm(y, mu1_cond, sqrt(sigma1_cond + 1 / x2), log = TRUE) + dlnorm(x2, mu2, sqrt(sigma2), log = TRUE))


      rbind(
        prob,
        mu1_cond_update * prob,
        log(x2) * prob,
        (sigma1_cond_update + mu1_cond_update**2) * prob,
        (mu1_cond_update * log(x2)) * prob,
        (mu1_cond_update * log(x2)) * prob,
        (log(x2)**2) * prob
      )
    }

    # val <- cubintegrate(f, c(exp(mu2-6*sqrt(sigma2))), c(exp(mu2+6*sqrt(sigma2))), fDim = 7, nVec = 1000)$integral
    val <- cubintegrate(f, c(0), c(Inf), fDim = 7, nVec = 1000)$integral
    ft <- matrix(val[2:3] / val[1], 2, 1)
    Qt <- matrix(val[4:7], 2, 2) / val[1] - ft %*% t(ft)
  } else {
    f0 <- ft
    S0 <- ginv(Qt)

    d1.log.like <- function(x) {
      mu <- x[1]
      tau <- exp(x[2])
      log.tau <- x[2]

      # -0.5*tau*((y-mu)**2)+0.5*log.tau
      c(
        tau * (y - mu),
        -0.5 * tau * ((y - mu)**2) + 0.5
      ) - S0 %*% (x - f0)
    }

    d2.inv <- function(x) {
      mu <- x[1]
      tau <- exp(x[2])
      log.tau <- x[2]

      cross <- tau * (y - mu)

      mat <- matrix(
        c(
          -tau,
          cross,
          cross,
          -0.5 * tau * ((y - mu)**2)
        ),
        2, 2
      ) - S0

      return(mat)
    }
    mode <- f_root(d1.log.like, d2.inv, f0)$root
    H <- d2.inv(mode)
    S <- ginv(-H)
    ft <- mode
    Qt <- S
  }

  return(list("ft" = ft, "Qt" = Qt))
}

#### Multi normal ####

#' convert_multi_NG_Normal
#'
#' Calculate the parameters of the Normal-Gamma that best approximates the given Multivariate Normal distribution. The distribution obtained for each outcome is marginal.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Normal to the Normal-Gamma.
#' In this approach, we suppose that the first entry of the multivariate normal represents the mean of the observed data and the second represent the log variance.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @return The parameters of the conjugated distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Normal outcome}
convert_multi_NG_Normal <- function(ft, Qt, parms) {
  k <- length(ft)
  r <- -3 / 2 + sqrt(9 / 4 + 2 * k)
  mu_index <- parms$mu_index
  var_index <- parms$var_index

  ft <- matrix(ft, k, 1)
  ft_mean <- ft[mu_index]
  ft_var <- ft[var_index]

  Qt_diag <- diag(Qt)
  Qt_mean <- Qt_diag[mu_index]
  Qt_var <- Qt_diag[var_index]

  if (r == 1) {
    Qt_cov <- Qt[1, 2]
  } else {
    Qt_cov <- diag(Qt[mu_index, var_index])
  }

  mu0 <- ft_mean + Qt_cov
  c0 <- exp(-ft_var - Qt_var / 2) / (Qt_mean + 1e-40)
  helper <- -3 + 3 * sqrt(1 + 2 * Qt_var / 3)
  # helper=Qt[2,2]
  alpha <- 1 / helper
  beta <- alpha * exp(-ft_var - Qt_var / 2)

  vec_par <- c(rbind(mu0, c0, alpha, beta))

  return(vec_par)
}

convert_multi_Normal_NG <- function(conj_distr, parms = list()) {
  r <- dim(conj_distr)[2] / 4
  k <- r + r * (r + 1) / 2

  mu0 <- conj_distr[, seq(1, r * 4 - 1, 4)]
  c0 <- conj_distr[, seq(1, r * 4 - 2, 4) + 1]
  alpha0 <- conj_distr[, seq(1, r * 4 - 3, 4) + 2]
  beta0 <- conj_distr[, seq(1, r * 4, 4) + 3]

  f1 <- mu0
  f2 <- digamma(alpha0) - log(beta0 + 1e-40)
  q1 <- beta0 / (c0 * alpha0)
  q2 <- trigamma(alpha0)
  q12 <- 0
  if (r == 1) {
    ft <- c(f1, f2)
    Qt <- matrix(c(q1, q12, q12, q2), byrow = F, ncol = 2)
  } else {
    ft <- c(f1, diag(f2)[upper.tri(diag(f2))])
    Qt <- diag(k)
    diag(Qt) <- c(q1, diag(q2)[upper.tri(diag(q2))])
  }
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_multi_NG_correl
#'
#'
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#' @importFrom Rfast lower_tri upper_tri lower_tri.assign upper_tri.assign
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Normal outcome}
#'
#' @details
#'
#' For evaluating the posterior parameters, we iterate over the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, updating one coordinate of the observation at a time.
#'
#' Since the original methodology requires a linear structure, a linearization is applied to the condinal mean and variance at each step. See \insertCite{ArtigoMultinormal;textual}{kDGLM}
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the detail about the methodology, see \insertCite{ArtigokParametrico;textual}{kDGLM} and \insertCite{ArtigoMultinormal;textual}{kDGLM}.
#'
#' @references
#'    \insertAllCited{}
update_multi_NG_correl <- function(conj_prior, ft, Qt, y, parms) {

  ######### TESTE INTERNO #########
  # lower.index=outcome$parms$lower.index
  # upper.index=outcome$parms$upper.index
  #
  # f=function(x){
  #
  #   mu <- x
  #   rho <- matrix(0, r, r)
  #   rho[upper.index] <- rho[lower.index] <- x[cor_index]
  #   diag(rho)=x[var_index]
  #   A=eigen(rho)
  #   vals=exp(A$values)
  #   Q=A$vectors
  #   Sigma <- Q %*% diag(vals) %*% t(Q)
  #   return(c(Sigma))
  # }
  # ft_up=rnorm(k)
  # f(ft_up)
  #
  # x=ft_up
  # mu <- x
  # rho[upper.index] <- rho[lower.index] <- x[cor_index]
  # diag(rho)=x[var_index]
  # A=eigen(rho)
  # vals=exp(A$values)
  # Q=A$vectors
  # Sigma <- Q %*% diag(vals) %*% t(Q)
  #
  # dSigma=(matrix(vals,r,r)-matrix(vals,r,r,byrow=TRUE))/
  #   ((matrix(A$values,r,r)-matrix(A$values,r,r,byrow=TRUE)))
  # diag(dSigma)=vals
  # B=matrix(0,r,r)
  # B[1,1]=1
  #
  #
  # D=matrix(0,r,r)
  # for(i in r){
  #   for(j in r){
  #     D=D+A$values[i]*Q[,i]%*%t(Q[j,])# %*% B %*% Q[j,]%*%t(Q[j,]) * dSigma[i,j]
  #   }
  # }
  # D
  # rho
  #
  # dS=t(calculus::derivative(f,var=ft_up))
  # matrix(dS[var_index[1],],r,r)
  # D
  #
  #
  # A=eigen(matrix(0,r,r))
  # vals=exp(A$values)
  # Q=A$vectors
  # Q %*% diag(vals) %*% t(Q)
  #################################

  mu_index <- parms$mu_index
  var_index <- parms$var_index
  cor_index <- parms$cor_index
  upper.index <- parms$upper.index
  lower.index <- parms$lower.index
  alt_method <- parms$alt_method

  r <- length(y)
  k <- r + r * (r + 1) / 2
  vec_r <- 1:(r**2)

  ft_up <- ft
  Qt_up <- Qt

  A <- matrix(0, k, 2)
  A[1, 1] <- A[r + 1, 2] <- 1
  ft_now <- ft_up[c(1, r + 1)]
  Qt_now <- Qt_up[c(1, r + 1), c(1, r + 1)]

  if (parms$alt_method) {
    post <- update_NG_alt(param, ft_now, Qt_now, y[1])
  } else {
    param <- convert_NG_Normal(ft_now, Qt_now)
    up_param <- update_NG(param, ft_now, Qt_now, y[1])
    post <- convert_Normal_NG(up_param)
  }
  print(Qt_now)

  ft_post <- post$ft
  Qt_post <- post$Qt

  At <- Qt_up[, c(1, r + 1)] %*% ginv(Qt_now)
  At_t <- t(At)
  ft_up <- ft_up + At %*% (ft_post - ft_now)
  Qt_up <- Qt_up + At %*% (Qt_post - Qt_now) %*% At_t

  if (r > 1) {
    for (i in 2:r) {
      {
      x <- c(ft_up)
      rho <- matrix(0, r, r)
      rho[upper.index] <- rho[lower.index] <- x[cor_index]
      var <- diag(exp(-x[var_index] / 2))
      p <- 1 / (1 + exp(-rho))
      rho <- 2 * p - 1
      diag(rho) <- 1
      Sigma <- var %*% rho %*% var
      # diag(rho)=diag(var)
      # Sigma=rho%*%t(rho)
      # Sigma=crossprod(transpose(rho))
      # print(eigen(Sigma))

      Sigma_rho <- Sigma[i, 1:(i - 1)]
      Sigma_part <- Sigma[1:(i - 1), 1:(i - 1)]

      if (all(Sigma_part == 0)) {
        ft_now <- x[c(i, i + r)]
        Qt_now <- Qt_up[c(i, i + r), c(i, i + r)]
      } else {
        S <- ginv(Sigma_part)
        e <- (y[1:(i - 1)] - x[1:(i - 1)])
        Sigma_S <- c(Sigma_rho %*% S)
        mu_bar <- x[i] + Sigma_S %*% e
        S_bar <- Sigma[i, i] - Sigma_S %*% Sigma_rho

        A <- matrix(0, k, 2)
        A[1:i, 1] <- c(-Sigma_S, 1)

        dx <- array(0, c(r, r, r * (r + 1) / 2))
        for (j in 1:i) {
          dx[j, j, j] <- -0.5 * var[j, j]
        }
        ref_p <- c(p)[upper.index]
        dx[vec_r[upper.index] + c(0:(k - 2 * r - 1)) * (r**2) + (r**3)] <-
          dx[vec_r[lower.index] + c(0:(k - 2 * r - 1)) * (r**2) + (r**3)] <-
          2 * ref_p * (1 - ref_p)

        dSigma <- array(NA, c(r, r, r * (r + 1) / 2))
        aux_1 <- rho %*% var

        dSigma[, , 1:r] <- array_mult_left(dx[, , 1:r], aux_1)

        dSigma[, , 1:r] <- dSigma[, , 1:r] + array_transp(dSigma[, , 1:r])
        dSigma[, , -(1:r)] <- dx[, , -(1:r), drop = FALSE] %>%
          array_mult_left(var) %>%
          array_mult_right(var)

        dSigma_part <- dSigma[(1:(i - 1)), (1:(i - 1)), , drop = FALSE]
        dSigma_rho <- dSigma[i, (1:(i - 1)), ] %>% matrix(i - 1, r * (r + 1) / 2)

        dSigma_p1 <- -array_mult_right(dSigma_part, S)
        dSigma_p1 <- array_mult_left(dSigma_p1, S)
        dSigma_p1 <- array_collapse_left(dSigma_p1, e)
        dSigma_p1 <- c(Sigma_rho %*% dSigma_p1)

        dSigma_p2 <- c(c(S %*% e) %*% dSigma_rho)
        A[-(1:r), 1] <- dSigma_p1 + dSigma_p2

        dSigma_p1 <- -array_collapse_right(dSigma_part, Sigma_S)
        dSigma_p1 <- c(Sigma_S %*% dSigma_p1)

        helper_p2 <- c(S %*% Sigma_rho)
        dSigma_p2 <- 2 * c(helper_p2 %*% dSigma_rho)
        A[-(1:r), 2] <- -(dSigma[i, i, ] - dSigma_p1 - dSigma_p2) / c(S_bar)

        ft_now <- c(mu_bar, -log(S_bar))
        Qt_now <- t(A) %*% Qt_up %*% A
      }
    }
    # {f=function(x){
    #
    #   mu <- x
    #   rho <- matrix(0, r, r)
    #   rho <- lower_tri.assign(rho, x[cor_index], diag = FALSE)
    #   rho <- upper_tri.assign(rho, x[cor_index], diag = FALSE)
    #   var <- diag(exp(-x[var_index] / 2))
    #   p <- 1 / (1 + exp(-rho))
    #   rho <- 2 * p - 1
    #   diag(rho)=1
    #   Sigma <- var %*% rho %*% var
    #
    #   S=ginv(Sigma[1:(i-1),1:(i-1)])
    #   mu_bar=mu[i]+Sigma[i,1:(i-1)]%*%S%*%(y[1:(i-1)]-mu[1:(i-1)])
    #   S_bar=Sigma[i,i]-Sigma[i,1:(i-1)]%*%S%*%Sigma[1:(i-1),i]
    #   # Sigma_S=solve(Sigma[1:(i-1),1:(i-1)],Sigma[1:(i-1),i])
    #   # mu_bar=mu[i]+Sigma_S%*%(y[1:(i-1)]-mu[1:(i-1)])
    #   # S_bar=Sigma[i,i]-Sigma_S%*%Sigma[1:(i-1),i]
    #   return(c(mu_bar,-log(S_bar)))
    # }
    #
    # A=t(calculus::derivative(f,var=ft_up))
    # ft_now=f(ft_up)
    # Qt_now=t(A)%*%Qt_up%*%A
    # }

    if (parms$alt_method) {
      post <- update_NG_alt(param, ft_now, Qt_now, y[i])
    } else {
      param <- convert_NG_Normal(ft_now, Qt_now)
      up_param <- update_NG(param, ft_now, Qt_now, y[i])
      post <- convert_Normal_NG(up_param)
    }
    print(Qt_now)

    ft_post <- post$ft
    Qt_post <- post$Qt

    At <- Qt_up %*% A %*% ginv(Qt_now)
    At_t <- t(At)

    ft_up <- ft_up + At %*% (ft_post - ft_now)
    Qt_up <- Qt_up + At %*% (Qt_post - Qt_now) %*% At_t    }
  }



  return(list("ft" = ft_up, "Qt" = Qt_up))
}

#' multi_normal_gamma_pred
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribution (Normal-Gamma) of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time.
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
#' @importFrom stats qt dt
#' @importFrom Rfast data.frame.to_matrix
#' @keywords internal
#' @family {auxiliary functions for a Normal outcome}
multi_normal_gamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  t <- if.null(dim(conj_param)[1], 1)
  r <- if.null(dim(conj_param)[2], length(conj_param)) / 4
  k <- r + r * (r + 1) / 2

  conj_param <- data.frame.to_matrix(conj_param)
  if (dim(conj_param)[2] == 1) {
    conj_param <- conj_param %>% t()
  }

  mu0 <- conj_param[, seq(1, r * 4 - 1, 4), drop = FALSE]
  c0 <- conj_param[, seq(1, r * 4 - 2, 4) + 1, drop = FALSE]
  alpha0 <- conj_param[, seq(1, r * 4 - 3, 4) + 2, drop = FALSE]
  beta0 <- conj_param[, seq(1, r * 4, 4) + 3, drop = FALSE]

  nu <- 2 * alpha0
  sigma2 <- (beta0 / alpha0) * (1 + 1 / c0)

  if (pred.flag) {
    pred <- mu0 %>% t()
    if(r==1){
      var.pred=array(sigma2,c(1,1,t))
    }else{
      var.pred <- array(0,c(r,r,t))
      for(i in 1:t){
        diag(var.pred[,,i])=sigma2[,i]
      }
    }

    icl.pred <- (qt((1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu0) %>% t()
    icu.pred <- (qt(1 - (1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu0) %>% t()
  } else {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }
  if (like.flag) {
    outcome <- matrix(outcome, t, r)
    log.like <- colSums(dt((outcome - mu0) * sqrt(alpha0 / beta0), nu, log = TRUE))
  } else {
    log.like <- NULL
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}
