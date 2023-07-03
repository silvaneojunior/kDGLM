#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G Array: A 3D-array representing the transition matrix of the model at each time.
#' @param G_labs Matrix: A character matrix containing the type associated with each value in G.
#'
#' @return List: The smoothed mean (mts) and covariance (Cts) of the latent variables at each time. Their dimension follows, respectively, the dimensions of mt and Ct.
#'
#' @importFrom Rfast transpose
#' @keywords internal
#'
#' @details
#'
#' For the models covered in this package, we always assume that the latent states have Multivariate Normal distribution. With that assumption, we can use Kalman Smoother algorithm to calculate the posterior of the states at each time, given everything that has been observed (assuming that we already know the filtered distribution of the states).
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about the algorithm implemented see \insertCite{ArtigokParametrico;textual}{kDGLM}, \insertCite{Petris-DLM;textual}{kDGLM}, chapter 2, \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 4, and \insertCite{Kalman_filter_origins;textual}{kDGLM}.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{analytic_filter}}
#' @references
#'    \insertAllCited{}
generic_smoother <- function(mt, Ct, at, Rt, G, G_labs) {
  T <- dim(mt)[2]
  mts <- mt
  Cts <- Ct
  if (T > 1) {
    for (t in (T - 1):1) {
      mt_now <- mt[, t, drop = FALSE]
      Ct_now <- Ct[, , t]
      G_step <- calc_current_G(mt_now, Ct_now, G[, , t + 1], G_labs)$G
      restricted_Rt <- Rt[, , t + 1]
      restricted_Ct <- Ct[, , t]

      simple_Rt_inv <- restricted_Ct %*% transpose(G_step) %*% ginv(restricted_Rt)
      simple_Rt_inv_t <- transpose(simple_Rt_inv)

      mts[, t] <- mt_now + simple_Rt_inv %*% (mts[, t + 1] - at[, t + 1])
      Cts[, , t] <- restricted_Ct + simple_Rt_inv %*% (Cts[, , t + 1] - restricted_Rt) %*% simple_Rt_inv_t
    }
  }
  return(list("mts" = mts, "Cts" = Cts))
}

#' analytic_filter
#'
#' Fit the model given the observed value and the model parameters.
#'
#' @param outcomes List: The observed data. It should contain objects of the class dlm_distr.
#' @param a1 Vector: The prior mean for the latent vector.
#' @param R1 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. Its dimension should be n x k x t, where n is the number of latent variables, k is the number of linear predictors in the model and t is the time series length.
#' @param FF_labs Matrix: A character matrix containing the type associated with each value in FF.
#' @param G Array: A 3D-array containing the state evolution matrix at each time.
#' @param G_labs Matrix: A character matrix containing the type associated with each value in G.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length.
#' @param h Matrix: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). Its dimension should be n x t, where t is the length of the series and n is the number of latent states.
#' @param H Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param p_monit numeric (optional): The prior probability of changes in the latent space variables that are not part of it's dynamic.
#' @param c_monit numeric (optional, if p_monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
#' @param monitoring Vector: A vector of flags indicating which latent variables should be monitored.
#'
#' @importFrom Rfast transpose
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x t.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x t.
#'    \item at Matrix: The one-step-ahead mean of the latent variables at each time. Dimensions are n x t.
#'    \item Rt Array: A 3D-array containing the one-step-ahead covariance matrix for latent variables at each time. Dimensions are n x n x t.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x t.
#'    \item Qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x m x t.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item G_labs Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values) when there is no monitoring. When monitoring for abnormalities, the value in times where abnormalities were detected is increased.
#'    \item h Array: The same as the argument (same values).
#'    \item H Array: The same as the argument (same values).
#'    \item W Array: A 3D-array containing the effective covariance matrix of the noise for each time, i.e., considering both H and D. Its dimension are the same as H and D.
#'    \item monitoring Vector: The same as the argument (same values).
#'    \item outcome List: The same as the argument outcome (same values).
#'    \item pred_names Vector: The names of the linear predictors.
#' }
#'
#' @keywords internal
#' @details
#'
#' For the models covered in this package, we always use the approach described in \insertCite{ArtigokParametrico;textual}{kDGLM}, including, in particular, the filtering algorithm presented in that work.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about the algorithm implemented see \insertCite{ArtigokParametrico;textual}{kDGLM}, \insertCite{Petris-DLM;textual}{kDGLM}, chapter 2, \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 4, and \insertCite{Kalman_filter_origins;textual}{kDGLM}.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{generic_smoother}}
#' @references
#'    \insertAllCited{}
analytic_filter <- function(outcomes, a1 = 0, R1 = 1, FF, FF_labs, G, G_labs, D, h, H, p_monit = NA, c_monit = 1, monitoring = FALSE) {
  # Defining quantities
  T <- dim(FF)[3]
  r <- sum(sapply(outcomes, function(x) {
    x$r
  }))
  n <- dim(FF)[1]
  k <- dim(FF)[2]

  pred_names <- colnames(FF)
  D_flags <- (D == 0) | array((t(monitoring) %>% crossprod(., .)) == 0, c(n, n, T))
  D <- ifelse(D == 0, 1, D)
  D_inv <- 1 / D
  a1 <- matrix(a1, n, 1)
  R1 <- R1
  mt <- matrix(NA, nrow = n, ncol = T)
  Ct <- array(NA, dim = c(n, n, T))
  Rt <- array(NA, dim = c(n, n, T))
  W <- array(NA, dim = c(n, n, T))
  ft <- matrix(NA, nrow = k, ncol = T)
  at <- matrix(NA, nrow = n, ncol = T)
  Qt <- array(NA, dim = c(k, k, T))
  monitoring <- array(monitoring, c(n))


  last_m <- a1
  last_C <- R1
  models <- list()
  D_mult <- list("null_model" = 1, "alt_model" = 100)
  H_add <- list("null_model" = 0, "alt_model" = 0.001)
  monit_win <- 1

  for (outcome_name in names(outcomes)) {
    if (!inherits(outcomes[[outcome_name]], "dlm_distr")) {
      stop(paste0("Error: Outcome contains is not of the right class Expected a dlm_distr object, got a ", class(outcomes[[outcome_name]]), " object."))
    }
    param_names <- outcomes[[outcome_name]]$param_names
    outcomes[[outcome_name]]$conj_prior_param <- matrix(NA, T, length(param_names)) %>% as.data.frame()
    names(outcomes[[outcome_name]]$conj_prior_param) <- param_names

    outcomes[[outcome_name]]$log.like.null <- rep(NA, T)
    outcomes[[outcome_name]]$log.like.alt <- rep(-Inf, T)
    outcomes[[outcome_name]]$alt.flags <- rep(0, T)


    pred_index <- match(outcomes[[outcome_name]]$pred_names, pred_names)
    k_i <- length(pred_index)

    outcomes[[outcome_name]]$ft <- matrix(NA, nrow = k_i, ncol = T)
    outcomes[[outcome_name]]$Qt <- array(NA, dim = c(k_i, k_i, T))
  }

  c <- c_monit
  p <- if.na(p_monit, 0)
  threshold <- log(c_monit) + log(p) - log(1 - p)
  last_C_D <- last_C

  for (t in 1:T) {
    model_list <- c("null_model")
    if (!is.na(p_monit)) {
      model_list <- c("alt_model", model_list)
    }
    for (model in model_list) {
      D_p <- D_inv[, , t]
      D_p[!D_flags[, , t]] <- D_p[!D_flags[, , t]] * D_mult[[model]]

      next_step <- one_step_evolve(last_m, last_C, G[, , t] %>% matrix(n, n), G_labs, D_p**0, h[, t], H[, , t] + diag(n) * H_add[[model]] + last_C_D * (D_p - 1))

      models[[model]] <- list(
        "at" = next_step$at, "Rt" = next_step$Rt,
        "at_step" = next_step$at, "Rt_step" = next_step$Rt,
        "W" = next_step$W
      )
    }
    FF_step <- matrix(FF[, , t], n, k, dimnames = list(NULL, pred_names))
    for (outcome_name in names(outcomes)) {
      outcome <- outcomes[[outcome_name]]
      pred_index <- match(outcome$pred_names, pred_names)
      # k_i <- length(pred_index)
      # FF_step <- matrix(FF[, pred_index, t], n, k_i)

      offset_step <- outcome$offset[t, ]
      na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)) | any(is.na(outcome$outcome[t, ])))
      for (model in model_list) {
        at_step <- models[[model]]$at_step
        Rt_step <- models[[model]]$Rt_step

        lin_pred <- calc_lin_pred(at_step, Rt_step, FF_step, FF_labs)

        ft_canom <- lin_pred$ft[pred_index, , drop = FALSE]
        Qt_canom <- lin_pred$Qt[pred_index, pred_index, drop = FALSE]
        if (outcome$convert_canom_flag) {
          ft_canom <- outcome$convert_mat_canom %*% ft_canom
          Qt_canom <- outcome$convert_mat_canom %*% Qt_canom %*% transpose(outcome$convert_mat_canom)
        }
        if (!na.flag) {
          offset_pred <- outcome$apply_offset(ft_canom, Qt_canom, offset_step)
          ft_canom <- offset_pred$ft
          Qt_canom <- offset_pred$Qt
        }

        conj_prior <- outcome$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)
        log.like <- NULL
        if (!is.na(p_monit)) {
          log.like <- outcome$calc_pred(conj_prior, outcome$outcome[t, ], parms = outcome$parms, pred_cred = NA)$log.like
        }

        models[[model]] <- list(
          "at" = models[[model]]$at, "Rt" = models[[model]]$Rt,
          "at_step" = at_step, "Rt_step" = Rt_step,
          "W" = models[[model]]$W,
          "ft_step" = ft_canom, "Qt_step" = Qt_canom,
          "FF_step" = lin_pred$FF[, pred_index, drop = FALSE],
          "conj_prior" = conj_prior,
          "log_like" = log.like
        )
      }
      model <- models$null_model
      if (!is.na(p_monit)) {
        outcomes[[outcome_name]]$log.like.null[t] <- models$null_model$log_like
        outcomes[[outcome_name]]$log.like.alt[t] <- models$alt_model$log_like
        bayes_factor <- sum(outcomes[[outcome_name]]$log.like.null[t:(t - monit_win + 1)] - outcomes[[outcome_name]]$log.like.alt[t:(t - monit_win + 1)]) %>% if.nan(0)

        if (monit_win > 0) {
          if (bayes_factor < threshold) {
            model <- models$alt_model
            conj_prior <- models$alt_model$conj_prior
            D_inv[, , t] <- D_inv[, , t] * D_mult$alt_model
            H[, , t] <- H[, , t] * H_add$alt_model
            monit_win <- -6
            outcomes[[outcome_name]]$alt.flags[t] <- 1
          } else if (bayes_factor > 0) {
            monit_win <- 0
          }
        }
        monit_win <- monit_win + 1
      }

      mt_step <- model$at_step
      Ct_step <- model$Rt_step
      ft_step <- model$ft_step
      Qt_step <- model$Qt_step
      norm_post <- list(ft = ft_step, Qt = Qt_step)
      if (!na.flag) {
        conj_prior <- outcome$conj_prior(ft_step, Qt_step, parms = outcome$parms)
        conj_post <- outcome$update(conj_prior, ft = ft_step, Qt = Qt_step, y = outcome$outcome[t, ], parms = outcome$parms)

        if (outcome$alt_method) {
          norm_post <- conj_post
        } else {
          norm_post <- outcome$conj_post(conj_post, parms = outcome$parms)
        }

        ft_star <- norm_post$ft <- norm_post$ft
        Qt_star <- norm_post$Qt <- norm_post$Qt
        error_ft <- (ft_star - ft_step)
        error_Qt <- (Qt_star - Qt_step)
        if (outcome$convert_canom_flag) {
          error_ft <- outcome$convert_mat_default %*% error_ft
          error_Qt <- outcome$convert_mat_default %*% error_Qt %*% transpose(outcome$convert_mat_default)
          Qt_step <- outcome$convert_mat_default %*% Qt_step %*% transpose(outcome$convert_mat_default)
        }
        At <- Ct_step %*% model$FF %*% ginv(Qt_step)
        models[["null_model"]]$at_step <- mt_step <- mt_step + At %*% error_ft
        models[["null_model"]]$Rt_step <- Ct_step <- Ct_step + At %*% error_Qt %*% t(At)
        last_C_D <- Ct_step
      }

      at[, t] <- model$at
      Rt[, , t] <- model$Rt

      lin_pred_ref <- calc_lin_pred(model$at, model$Rt, FF_step, FF_labs)
      ft[, t] <- lin_pred_ref$ft
      Qt[, , t] <- lin_pred_ref$Qt

      for (outcome_name in names(outcomes)) {
        outcome <- outcomes[[outcome_name]]
        pred_index <- match(outcome$pred_names, pred_names)

        offset_step <- outcome$offset[t, ]
        na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)))

        ft_canom <- lin_pred_ref$ft[pred_index, , drop = FALSE]
        Qt_canom <- lin_pred_ref$Qt[pred_index, pred_index, drop = FALSE]

        lin_pred <- list(ft = ft_canom, Qt = Qt_canom)

        if (outcome$convert_canom_flag) {
          ft_canom <- outcome$convert_mat_canom %*% ft_canom
          Qt_canom <- outcome$convert_mat_canom %*% Qt_canom %*% transpose(outcome$convert_mat_canom)
        }
        if (!na.flag) {
          offset_pred <- outcome$apply_offset(ft_canom, Qt_canom, offset_step)
          ft_canom <- offset_pred$ft
          Qt_canom <- offset_pred$Qt
        }
        conj_prior <- outcome$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)
        outcomes[[outcome_name]]$conj_prior_param[t, ] <- conj_prior
      }
      mt[, t] <- last_m <- mt_step
      Ct[, , t] <- last_C <- Ct_step
      W[, , t] <- model$W
    }
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "FF" = FF, "FF_labs" = FF_labs,
    "G" = G, "G_labs" = G_labs,
    "D" = D, "h" = h, "H" = H, "W" = W,
    "monitoring" = monitoring,
    "outcomes" = outcomes, "pred_names" = pred_names
  )
  return(result)
}

calc_current_G <- function(m0, C0, G, G_labs) {
  n <- if.null(dim(G)[1], 1)
  G_now <- G %>% as.matrix()

  drift <- m0 * 0
  flag_na <- rowSums(is.na(G_now)) >= 1
  index_na <- (1:n)[flag_na]
  if (any(flag_na)) {
    for (index_row in index_na) {
      index_col <- (1:n)[is.na(G_now[index_row, ])]
      flags_const <- G_labs[index_row, index_col] == "constrained"
      if (any(flags_const)) {
        index_const <- index_col[flags_const]
        rho <- tanh(m0[index_const + 1])
        G_now[index_row, index_const] <- rho
        G_now[index_row, index_const + 1] <- (1 + rho) * (1 - rho) * m0[index_const]
        drift[index_row] <- drift[index_row] - sum((1 + rho) * (1 - rho) * m0[index_col] * m0[index_col + 1])
      }

      flags_free <- G_labs[index_row, index_col] == "free"
      if (any(flags_free)) {
        index_free <- index_col[flags_free]
        G_now[index_row, index_free] <- m0[index_free + 1]
        G_now[index_row, index_free + 1] <- m0[index_free]
        drift[index_row] <- drift[index_row] - sum(m0[index_free] * m0[index_free + 1])
      }


      m1 <- m0
      C1 <- C0
      flags_kl <- G_labs[index_row, index_col] == "kl"
      if (any(flags_kl)) {
        index_kl <- index_col[flags_kl]
        for (i in index_kl) {
          m1[i] <- m0[i] * m0[i + 1] + C0[i, i + 1]
          C1[i, i] <- C0[i, i] * C0[i + 1, i + 1] +
            C0[i, i + 1]**2 +
            2 * m0[i] * m0[i + 1] * C0[i, i + 1] +
            (m0[i]**2) * C0[i + 1, i + 1] + (m0[i + 1]**2) * C0[i, i]
          cov_vec <- c(
            C0[i, i + 1] * m0[i] + C0[i, i] * m0[i + 1],
            C0[i, i + 1] * m0[i + 1] + C0[i + 1, i + 1] * m0[i]
          )
          C1[-i, i] <-
            C1[i, -i] <-
            cov_vec %*% ginv(C0[i:(i + 1), i:(i + 1)]) %*% C0[i:(i + 1), -i]
          G_now[i, i] <- 1
        }
        C1 <- G_now %*% C1 %*% transpose(G_now)
        G_now <- create_G(C0, C1)
        drift <- drift - G_now %*% m0 + m1

        # n_kl=sum(flags_kl)
        # m1[index_row]=0
        # C1[index_row, index_row]=0
        # C1[-index_row, index_row]=0
        # C1[index_row, -index_row]=0
        #
        # cross_prod=m0%*%t(m0)+C0
        #
        # m1[index_row] <- sum(diag(cross_prod[index_kl,index_kl+1]))
        # C1[index_row, index_row] <- sum(
        #   cross_prod[flags_kl,flags_kl]*cross_prod[flags_kl+1,flags_kl+1]+
        #     cross_prod[flags_kl,flags_kl+1]*cross_prod[flags_kl+1,flags_kl]
        # )
        #
        # if(n_kl>1){
        #   out_kl=index_kl[index_kl!=i]
        #
        #
        #   cross_prod[i,out_kl]*cross_prod[i+1,out_kl+1]+
        #     cross_prod[i,out_kl+1]*cross_prod[i+1,out_kl]
        #
        #   C1[index_row, index_row] <- C1[index_row, index_row]+
        #     C0[i, i] * C0[i + 1, i + 1] +
        #     C0[i, i + 1]**2 +
        #     2 * m0[i] * m0[i + 1] * C0[i, i + 1] +
        #     (m0[i]**2) * C0[i + 1, i + 1] + (m0[i + 1]**2) * C0[i, i]
        # }
      }
    }
  }

  m1 <- G_now %*% m0 + drift
  C1 <- G_now %*% C0 %*% transpose(G_now)

  list("m1" = m1, "C1" = C1, "G" = G_now, "drift" = drift)
}

one_step_evolve <- function(m0, C0, G, G_labs, D_inv, h, H) {
  G_vals <- calc_current_G(m0, C0, G, G_labs)

  # print('####################')
  # print(cbind(m0,G_vals$m1))
  # print(cbind(C0,G_vals$C1))
  G_now <- G_vals$G
  drift <- G_vals$drift
  at <- G_vals$m1 + h
  Pt <- G_vals$C1
  Rt <- as.matrix(D_inv * Pt) + H

  list("at" = at, "Rt" = Rt, "G" = G_now, "h" = h + drift, "W" = as.matrix((D_inv - 1) * Pt) + H)
}

calc_current_F <- function(at, Rt, FF, FF_labs) {
  pred.names <- colnames(FF)
  n <- dim(FF)[1]
  k <- dim(FF)[2]


  charge <- matrix(0, k, 1)
  at_mod <- at[, 1]
  count.na <- sum(is.na(FF))
  while (count.na > 0) {
    flag.na <- colSums(is.na(FF)) > 0
    index.na <- (1:k)[flag.na]
    for (index.pred in index.na) {
      flag.var <- is.na(FF[, index.pred])
      index.var <- (1:n)[flag.var]
      for (index.effect in index.var) {
        effect.name <- FF_labs[index.effect, index.pred]
        effect.vals <- FF[, effect.name == pred.names]
        if (any(is.na(effect.vals))) {
          break
        }
        FF[index.var, index.pred] <- sum(effect.vals * at_mod)
        FF[, index.pred] <- FF[, index.pred] + effect.vals * at_mod[index.var]
        charge[index.pred, 1] <- charge[index.pred, 1] - at_mod[index.var] * sum(effect.vals * at_mod)
      }
    }
    new.count.na <- sum(is.na(FF))
    if (count.na == new.count.na) {
      stop("Error: A circularity was detected in the specification of FF. Revise the definition of the linear predictors.")
    }
    count.na <- new.count.na
  }

  list("FF" = FF, "FF_diff" = charge)
}

calc_lin_pred <- function(at, Rt, FF, FF_labs) {
  FF_vals <- calc_current_F(at, Rt, FF, FF_labs)
  FF <- FF_vals$FF
  FF_diff <- FF_vals$FF_diff
  # print('#########################')
  # print(FF)
  # print(at)
  ft <- crossprod(FF, at) + FF_diff
  Qt <- as.matrix(crossprod(FF, Rt) %*% FF)
  list("ft" = ft, "Qt" = Qt, "FF" = FF, "FF_diff" = FF_diff)
}

format_ft <- function(ft, Qt, parms) {
  return(do.call(c, list(ft, Qt)))
}

format_param <- function(conj_param, parms) {
  if (is.null(dim(conj_param))) {
    r_star <- length(conj_param)
    r <- (-1 + sqrt(1 + 4 * r_star)) / 2
    t <- 1
    ft <- conj_param[1:r] %>% matrix(r, t)
    Qt <- conj_param[(r + 1):(r * r + r)] %>% array(c(r, r, t))
  } else {
    r_star <- dim(conj_param)[2]
    t <- dim(conj_param)[1]
    r <- (-1 + sqrt(1 + 4 * r_star)) / 2
    ft <- conj_param[, 1:r] %>%
      data.frame.to_matrix() %>%
      t()
    Qt <- conj_param[, (r + 1):(r * r + r)] %>%
      data.frame.to_matrix() %>%
      t() %>%
      array(c(r, r, t))
  }
  if (t == 1) {
    Qt <- matrix(Qt, r, r)
  }
  return(list("ft" = ft, "Qt" = Qt))
}

generic_param_names <- function(k) {
  c(
    paste("ft_", 1:k, sep = ""),
    paste("Qt_",
      c(matrix(1:k, k, k)),
      c(matrix(1:k, k, k, byrow = TRUE)),
      sep = ""
    )
  )
}
