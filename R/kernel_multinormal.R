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
  # parms=outcome$parms
  # ft_up=level$a1
  # Qt_up=level$R1
  # ft_star=array(NA,c(5,5))
  # Qt_star=array(NA,c(5,5,5))
  #
  # for(index in 1:5){
  #   ft=ft_up
  #   Qt=Qt_up
  # y=outcome$outcome[index,]

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

  ft_post <- post$ft
  Qt_post <- post$Qt
  # print('\n#################### first #######################')
  # print(Qt_now)
  # print(Qt_post)

  At <- Qt_up[, c(1, r + 1)] %*% ginv(Qt_now)
  At_t <- t(At)
  ft_up <- ft_up + At %*% (ft_post - ft_now)
  Qt_up <- Qt_up + At %*% (Qt_post - Qt_now) %*% At_t

  if (r > 1) {
    for (i in 2:r) {
      {      x <- c(ft_up)
        rho <- matrix(0, r, r)
        rho[upper.index] <- rho[lower.index] <- x[cor_index]
        sd <- diag(exp(-x[var_index] / 2))
        rho <- tanh(rho)
        diag(rho) <- 1
        Sigma <- sd %*% rho %*% sd
        # diag(rho)=diag(sd)
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
            dx[j, j, j] <- -0.5 * sd[j, j]
          }
          ref_rho <- c(rho)[upper.index]
          dx[vec_r[upper.index] + c(0:(k - 2 * r - 1)) * (r**2) + (r**3)] <-
            dx[vec_r[lower.index] + c(0:(k - 2 * r - 1)) * (r**2) + (r**3)] <-
            (1 + ref_rho) * (1 - ref_rho)

          dSigma <- array(NA, c(r, r, r * (r + 1) / 2))
          aux_1 <- rho %*% sd

          dSigma[, , 1:r] <- array_mult_left(dx[, , 1:r], aux_1)

          dSigma[, , 1:r] <- dSigma[, , 1:r] + array_transp(dSigma[, , 1:r])
          dSigma[, , -(1:r)] <- dx[, , -(1:r), drop = FALSE] |>
            array_mult_left(sd) |>
            array_mult_right(sd)

          dSigma_part <- dSigma[(1:(i - 1)), (1:(i - 1)), , drop = FALSE]
          dSigma_rho <- dSigma[i, (1:(i - 1)), ] |> matrix(i - 1, r * (r + 1) / 2)

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
      ###################################
      # {f=function(x){
      #
      #   mu <- x
      #   rho <- matrix(0, r, r)
      #   rho <- lower_tri.assign(rho, x[cor_index], diag = FALSE)
      #   rho <- upper_tri.assign(rho, x[cor_index], diag = FALSE)
      #   sd <- diag(exp(-x[var_index] / 2))
      #   rho <- tanh(rho)
      #   diag(rho)=1
      #   Sigma <- sd %*% rho %*% sd
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
      # A_test=t(calculus::derivative(f,sd=ft_up))
      # # ft_now=f(ft_up)
      # # Qt_now=t(A)%*%Qt_up%*%A
      # print(max(abs(A-A_test)))
      # }
      ####################################

      if (parms$alt_method) {
        post <- update_NG_alt(param, ft_now, Qt_now, y[i])
      } else {
        param <- convert_NG_Normal(ft_now, Qt_now)
        up_param <- update_NG(param, ft_now, Qt_now, y[i])
        post <- convert_Normal_NG(up_param)
      }

      ft_post <- post$ft
      Qt_post <- post$Qt
      # print('\n#################### second #######################')
      # print(Qt_now)
      # print(Qt_post)
      # print(Qt_up)

      At <- Qt_up %*% A %*% ginv(Qt_now)
      At_t <- t(At)

      ft_up <- ft_up + At %*% (ft_post - ft_now)
      Qt_up <- Qt_up + At %*% (Qt_post - Qt_now) %*% At_t
    }
    # ft_star[,index]=ft_up
    # Qt_star[,,index]=Qt_up
  }
  # }



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
    conj_param <- conj_param |> t()
  }

  mu0 <- conj_param[, seq(1, r * 4 - 1, 4), drop = FALSE]
  c0 <- conj_param[, seq(1, r * 4 - 2, 4) + 1, drop = FALSE]
  alpha0 <- conj_param[, seq(1, r * 4 - 3, 4) + 2, drop = FALSE]
  beta0 <- conj_param[, seq(1, r * 4, 4) + 3, drop = FALSE]

  nu <- 2 * alpha0
  sigma2 <- (beta0 / alpha0) * (1 + 1 / c0)

  if (pred.flag) {
    pred <- mu0 |> t()
    if (r == 1) {
      var.pred <- array(sigma2, c(1, 1, t))
    } else {
      var.pred <- array(0, c(r, r, t))
      for (i in 1:r) {
        var.pred[i, i, ] <- sigma2[, i]
      }
    }

    icl.pred <- (qt((1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu0) |> t()
    icu.pred <- (qt(1 - (1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu0) |> t()
  } else {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }
  if (like.flag) {
    outcome <- matrix(outcome, t, r)
    log.like <- colSums(dt((outcome - mu0) * sqrt(c0 * alpha0 / beta0), nu, log = TRUE) + log(c0 * alpha0 / beta0))
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
