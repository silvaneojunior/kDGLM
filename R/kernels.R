#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G Array: A 3D-array representing the transition matrix of the model at each time.
#' @param G.labs Matrix: A character matrix containing the type associated with each value in G.
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
generic_smoother <- function(mt, Ct, at, Rt, G, G.labs) {
  T_len <- dim(mt)[2]
  mts <- mt
  Cts <- Ct
  if (T_len > 1) {
    for (t in (T_len - 1):1) {
      mt.step <- mt[, t, drop = FALSE]
      Ct.step <- Ct[, , t]
      Rt.step <- Rt[, , t + 1]
      G.step <- calc_current_G(mt.step, Ct.step, G[, , t + 1], G.labs)$G
      # G.step <- G[, , t + 1]

      simple.Rt.inv <- Ct.step %*% transpose(G.step) %*% ginv(Rt.step)
      simple.Rt.inv.t <- transpose(simple.Rt.inv)

      mts[, t] <- mt.step + simple.Rt.inv %*% (mts[, t + 1] - at[, t + 1])
      Cts[, , t] <- Ct.step + simple.Rt.inv %*% (Cts[, , t + 1] - Rt.step) %*% simple.Rt.inv.t
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
#' @param FF.labs Matrix: A character matrix containing the type associated with each value in FF.
#' @param G Array: A 3D-array containing the state evolution matrix at each time.
#' @param G.labs Matrix: A character matrix containing the type associated with each value in G.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length.
#' @param h Matrix: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). Its dimension should be n x t, where t is the length of the series and n is the number of latent states.
#' @param H Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param p.monit numeric (optional): The prior probability of changes in the latent space variables that are not part of its dynamic.
#' @param c.monit numeric (optional, if p.monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting abnormalities.
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
#'    \item G.labs Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values) when there is no monitoring. When monitoring for abnormalities, the value in times where abnormalities were detected is increased.
#'    \item h Array: The same as the argument (same values).
#'    \item H Array: The same as the argument (same values).
#'    \item W Array: A 3D-array containing the effective covariance matrix of the noise for each time, i.e., considering both H and D. Its dimension are the same as H and D.
#'    \item monitoring Vector: The same as the argument (same values).
#'    \item outcomes List: The same as the argument outcomes (same values).
#'    \item pred.names Vector: The names of the linear predictors.
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
analytic_filter <- function(outcomes,
                            a1 = 0, R1 = 1,
                            FF, FF.labs,
                            G, G.labs,
                            D, h, H,
                            p.monit = NA, c.monit = 1,
                            monitoring = FALSE) {
  # Defining quantities
  T_len <- dim(FF)[3]
  n <- dim(FF)[1]
  k <- dim(FF)[2]

  pred.names <- colnames(FF)
  D.flags <- (D == 0) | array(crossprod(t(monitoring), t(monitoring)) == 0, c(n, n, T_len))
  D <- ifelse(D == 0, 1, D)
  D.inv <- 1 / D
  a1 <- matrix(a1, n, 1)
  R1 <- R1
  mt <- matrix(NA, nrow = n, ncol = T_len)
  Ct <- array(NA, dim = c(n, n, T_len))
  Rt <- array(NA, dim = c(n, n, T_len))
  W <- array(NA, dim = c(n, n, T_len))
  ft <- matrix(NA, nrow = k, ncol = T_len)
  at <- matrix(NA, nrow = n, ncol = T_len)
  Qt <- array(NA, dim = c(k, k, T_len))
  monitoring <- array(monitoring, c(n))
  null.rows.flag=rowSums(G.labs=='noise.disc',na.rm=TRUE)>0


  last.m <- a1
  last.C <- R1
  models <- list()
  D.mult <- list("null.model" = 1, "alt.model" = 100)
  H.add <- list("null.model" = diag(n) * 0, "alt.model" = diag(n) * 0.001)
  monit.win <- 1

  for (outcome.name in names(outcomes)) {
    if (!inherits(outcomes[[outcome.name]], "dlm_distr")) {
      stop(paste0("Error: Outcome contains is not of the right class Expected a dlm_distr object, got a ", class(outcomes[[outcome.name]]), " object."))
    }
    param.names <- outcomes[[outcome.name]]$param.names
    # outcomes[[outcome.name]]$conj.prior.param <- matrix(NA, T_len, length(param.names)) |> as.data.frame()
    # names(outcomes[[outcome.name]]$conj.prior.param) <- param.names

    outcomes[[outcome.name]]$log.like.null <- rep(NA, T_len)
    outcomes[[outcome.name]]$log.like.alt <- rep(-Inf, T_len)
    outcomes[[outcome.name]]$alt.flags <- rep(0, T_len)


    pred.index <- match(outcomes[[outcome.name]]$pred.names, pred.names)
    k_i <- length(pred.index)

    outcomes[[outcome.name]]$ft <- matrix(NA, nrow = k_i, ncol = T_len)
    outcomes[[outcome.name]]$Qt <- array(NA, dim = c(k_i, k_i, T_len))
  }

  c <- c.monit
  p <- if.na(p.monit, 0)
  threshold <- log(c.monit) + log(p) - log(1 - p)

  H.now=R1
  for (t in seq_len(T_len)) {
    model.list <- c("null.model")
    if (!is.na(p.monit)) {
      model.list <- c("alt.model", model.list)
    }
    for (model in model.list) {
      D.p.inv <- D.inv[, , t]
      D.p <- D[, , t]
      D.p.inv[!D.flags[, , t]] <- D.p.inv[!D.flags[, , t]] * D.mult[[model]]

      H.prev=H.now
      G.now=G[, , t] |> matrix(n, n)
      H.now=H[, , t] |> matrix(n, n)

      if(any(null.rows.flag)){
        weight=((t-1)/t)*diag(D.p)
        D.p.inv[null.rows.flag,null.rows.flag]=1
        diag(H.now)[null.rows.flag]=((weight*diag(H.prev))+
                                         ((1-weight)*(diag(last.C)+last.m**2)))[null.rows.flag]
      }

      next.step <- one_step_evolve(last.m, last.C,
                                   G.now, G.labs,
                                   D.p.inv, h[, t], H.now + H.add[[model]])


      models[[model]] <- list(
        "at" = next.step$at, "Rt" = next.step$Rt,
        "at.step" = next.step$at, "Rt.step" = next.step$Rt,
        "h" = next.step$h,
        "W" = next.step$W,
        "G" = next.step$G
      )
    }
    FF.step <- matrix(FF[, , t], n, k, dimnames = list(NULL, pred.names))
    for (outcome.name in names(outcomes)) {
      outcome <- outcomes[[outcome.name]]
      pred.index <- match(outcome$pred.names, pred.names)
      # k_i <- length(pred.index)
      # FF.step <- matrix(FF[, pred.index, t], n, k_i)

      offset.step <- outcome$offset[t, ]
      na.flag <- any(is.null(offset.step) | any(offset.step == 0) | any(is.na(offset.step)) | any(is.na(outcome$data[t, ])))
      for (model in model.list) {
        at.step <- models[[model]]$at.step
        Rt.step <- models[[model]]$Rt.step

        lin.pred <- calc_lin_pred(at.step, Rt.step, FF.step, FF.labs, pred.names)

        ft.canom <- lin.pred$ft[pred.index, , drop = FALSE]
        Qt.canom <- lin.pred$Qt[pred.index, pred.index, drop = FALSE]
        if (outcome$convert.canom.flag) {
          ft.canom <- outcome$convert.mat.canom %*% ft.canom
          Qt.canom <- outcome$convert.mat.canom %*% Qt.canom %*% transpose(outcome$convert.mat.canom)
        }
        if (!na.flag) {
          offset.pred <- outcome$apply_offset(ft.canom, Qt.canom, offset.step)
          ft.canom <- offset.pred$ft
          Qt.canom <- offset.pred$Qt
        }

        conj.prior <- outcome$conj_distr(ft.canom, Qt.canom, parms = outcome$parms)
        log.like <- NULL
        if (!is.na(p.monit)) {
          log.like <- outcome$calc_pred(conj.prior, outcome$data[t, ], parms = outcome$parms, pred.cred = NA)$log.like
        }

        models[[model]] <- list(
          "at" = models[[model]]$at, "Rt" = models[[model]]$Rt,
          "at.step" = at.step, "Rt.step" = Rt.step,
          "h" = models[[model]]$h,
          "W" = models[[model]]$W,
          "G" = models[[model]]$G,
          "ft.step" = ft.canom, "Qt.step" = Qt.canom,
          "FF.step" = lin.pred$FF[, pred.index, drop = FALSE],
          "conj.prior" = conj.prior,
          "log.like" = log.like
        )
      }
      model <- models$null.model
      if (!is.na(p.monit)) {
        outcomes[[outcome.name]]$log.like.null[t] <- models$null.model$log.like
        outcomes[[outcome.name]]$log.like.alt[t] <- models$alt.model$log.like
        bayes.factor <- sum(outcomes[[outcome.name]]$log.like.null[t:(t - monit.win + 1)] +
                              -outcomes[[outcome.name]]$log.like.alt[t:(t - monit.win + 1)]) |>
          if.nan(0)

        if (monit.win > 0) {
          if (bayes.factor < threshold) {
            model <- models$alt.model
            conj.prior <- models$alt.model$conj.prior
            D.inv[, , t] <- D.inv[, , t] * D.mult$alt.model
            H[, , t] <- H[, , t] * H.add$alt.model
            monit.win <- -6
            outcomes[[outcome.name]]$alt.flags[t] <- 1
          } else if (bayes.factor > 0) {
            monit.win <- 0
          }
        }
        monit.win <- monit.win + 1
      }

      mt.step <- model$at.step
      Ct.step <- model$Rt.step
      ft.step <- model$ft.step
      Qt.step <- model$Qt.step
      norm.post <- list(ft = ft.step, Qt = Qt.step)
      if (!na.flag) {
        conj.prior <- outcome$conj_distr(ft.step, Qt.step, parms = outcome$parms)
        conj.post <- outcome$update(conj.prior, ft = ft.step, Qt = Qt.step, y = outcome$data[t, ], parms = outcome$parms)

        if (outcome$alt.method) {
          norm.post <- conj.post
        } else {
          norm.post <- outcome$norm_distr(conj.post, parms = outcome$parms)
        }

        ft.star <- norm.post$ft <- norm.post$ft
        Qt.star <- norm.post$Qt <- norm.post$Qt
        error.ft <- (ft.star - ft.step)
        error.Qt <- (Qt.star - Qt.step)
        if (outcome$convert.canom.flag) {
          error.ft <- outcome$convert.mat.default %*% error.ft
          error.Qt <- outcome$convert.mat.default %*% error.Qt %*% transpose(outcome$convert.mat.default)
          Qt.step <- outcome$convert.mat.default %*% Qt.step %*% transpose(outcome$convert.mat.default)
        }
        At <- Ct.step %*% model$FF %*% ginv(Qt.step)
        models[["null.model"]]$at.step <- mt.step <- mt.step + At %*% error.ft
        models[["null.model"]]$Rt.step <- Ct.step <- Ct.step + At %*% error.Qt %*% t(At)
      }

      at[, t] <- model$at
      Rt[, , t] <- model$Rt

      lin.pred.ref <- calc_lin_pred(model$at, model$Rt, FF.step, FF.labs, pred.names)
      ft[, t] <- lin.pred.ref$ft
      Qt[, , t] <- lin.pred.ref$Qt

      for (outcome.name in names(outcomes)) {
        outcome <- outcomes[[outcome.name]]
        pred.index <- match(outcome$pred.names, pred.names)

        offset.step <- outcome$offset[t, ]
        na.flag <- any(is.null(offset.step) | any(offset.step == 0) | any(is.na(offset.step)))

        ft.canom <- lin.pred.ref$ft[pred.index, , drop = FALSE]
        Qt.canom <- lin.pred.ref$Qt[pred.index, pred.index, drop = FALSE]

        lin.pred <- list(ft = ft.canom, Qt = Qt.canom)

        if (outcome$convert.canom.flag) {
          ft.canom <- outcome$convert.mat.canom %*% ft.canom
          Qt.canom <- outcome$convert.mat.canom %*% Qt.canom %*% transpose(outcome$convert.mat.canom)
        }
        if (!na.flag) {
          offset.pred <- outcome$apply_offset(ft.canom, Qt.canom, offset.step)
          ft.canom <- offset.pred$ft
          Qt.canom <- offset.pred$Qt
        }
        conj.prior <- outcome$conj_distr(ft.canom, Qt.canom, parms = outcome$parms)
        # outcomes[[outcome.name]]$conj.prior.param[t, ] <- conj.prior
      }
      mt[, t] <- last.m <- mt.step
      Ct[, , t] <- last.C <- Ct.step
      # h[, t] <- model$h
      W[, , t] <- model$W
      # G[, , t] <- model$G
    }
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "FF" = FF, "FF.labs" = FF.labs,
    "G" = G, "G.labs" = G.labs,
    "D" = D, "h" = h, "H" = H, "W" = W,
    "monitoring" = monitoring, smooth = FALSE,
    "outcomes" = outcomes, "pred.names" = pred.names
  )
  return(result)
}

calc_current_G <- function(m0, C0, G, G.labs) {
  n <- length(m0)
  G.now <- G |> matrix(n,n)

  drift <- m0 * 0
  flag.na <- rowSums(is.na(G.now)) >= 1
  index.na <- seq_len(n)[flag.na]
  if (any(flag.na)) {
    for (index.row in index.na) {
      index.col <- seq_len(n)[is.na(G.now[index.row, ])]
      flags.const <- G.labs[index.row, index.col] == "constrained"
      if (any(flags.const)) {
        index.const <- index.col[flags.const]
        rho <- tanh(m0[index.const + 1])
        G.now[index.row, index.const] <- rho
        G.now[index.row, index.const + 1] <- (1 + rho) * (1 - rho) * m0[index.const]
        drift[index.row] <- drift[index.row] - sum((1 + rho) * (1 - rho) * m0[index.col] * m0[index.col + 1])
      }

      flags.free <- G.labs[index.row, index.col] == "free"
      if (any(flags.free)) {
        index.free <- index.col[flags.free]
        G.now[index.row, index.free] <- m0[index.free + 1]
        G.now[index.row, index.free + 1] <- m0[index.free]
        drift[index.row] <- drift[index.row] - sum(m0[index.free] * m0[index.free + 1])
      }


      m1 <- m0
      C1 <- C0
      flags.kl <- G.labs[index.row, index.col] == "kl"
      if (any(flags.kl)) {
        index.kl <- index.col[flags.kl]
        for (i in index.kl) {
          m1[i] <- m0[i] * m0[i + 1] + C0[i, i + 1]
          C1[i, i] <- C0[i, i] * C0[i + 1, i + 1] +
            C0[i, i + 1]**2 +
            2 * m0[i] * m0[i + 1] * C0[i, i + 1] +
            (m0[i]**2) * C0[i + 1, i + 1] + (m0[i + 1]**2) * C0[i, i]
          cov.vec <- c(
            C0[i, i + 1] * m0[i] + C0[i, i] * m0[i + 1],
            C0[i, i + 1] * m0[i + 1] + C0[i + 1, i + 1] * m0[i]
          )
          C1[-i, i] <-
            C1[i, -i] <-
            cov.vec %*% ginv(C0[i:(i + 1), i:(i + 1)]) %*% C0[i:(i + 1), -i]
          G.now[i, i] <- 1
        }
        C1 <- G.now %*% C1 %*% transpose(G.now)
        G.now <- create_G(C0, C1)
        drift <- drift - G.now %*% m0 + m1

        # n_kl=sum(flags.kl)
        # m1[index.row]=0
        # C1[index.row, index.row]=0
        # C1[-index.row, index.row]=0
        # C1[index.row, -index.row]=0
        #
        # cross_prod=m0%*%t(m0)+C0
        #
        # m1[index.row] <- sum(diag(cross_prod[index.kl,index.kl+1]))
        # C1[index.row, index.row] <- sum(
        #   cross_prod[flags.kl,flags.kl]*cross_prod[flags.kl+1,flags.kl+1]+
        #     cross_prod[flags.kl,flags.kl+1]*cross_prod[flags.kl+1,flags.kl]
        # )
        #
        # if(n_kl>1){
        #   out_kl=index.kl[index.kl!=i]
        #
        #
        #   cross_prod[i,out_kl]*cross_prod[i+1,out_kl+1]+
        #     cross_prod[i,out_kl+1]*cross_prod[i+1,out_kl]
        #
        #   C1[index.row, index.row] <- C1[index.row, index.row]+
        #     C0[i, i] * C0[i + 1, i + 1] +
        #     C0[i, i + 1]**2 +
        #     2 * m0[i] * m0[i + 1] * C0[i, i + 1] +
        #     (m0[i]**2) * C0[i + 1, i + 1] + (m0[i + 1]**2) * C0[i, i]
        # }
      }
    }
  }

  m1 <- G.now %*% m0 + drift
  C1 <- G.now %*% C0 %*% transpose(G.now)

  list("m1" = m1, "C1" = C1, "G" = G.now, "drift" = drift)
}

one_step_evolve <- function(m0, C0, G, G.labs, D.inv, h, H) {
  G.vals <- calc_current_G(m0, C0, G, G.labs)

  # print('####################')
  # print(cbind(m0,G.vals$m1))
  # print(cbind(C0,G.vals$C1))
  G.now <- G.vals$G
  drift <- G.vals$drift
  at <- G.vals$m1 + h
  Pt <- G.vals$C1
  W <- ((D.inv - 1) * Pt) + H
  Rt <- W+Pt

  list("at" = at, "Rt" = Rt, "G" = G.now, "h" = h + drift, "W" = W)
}

calc_current_F <- function(at, Rt, FF, FF.labs,pred.names) {
  n <- dim(FF)[1]
  k <- dim(FF)[2]

  charge <- matrix(0, k, 1)
  at.mod <- c(at)
  count.na <- sum(is.na(FF))
  if (any(is.na(FF) & FF.labs == "const")) {
    stop("Error: Unexpected NA values in the FF matrix.")
  }
  while (count.na > 0) {
    flag.na <- colSums(is.na(FF)) > 0
    index.na <- seq_len(k)[flag.na]
    for (index.pred in index.na) {
      flag.var <- is.na(FF[, index.pred])
      index.var <- seq_len(n)[flag.var]
      for (index.effect in index.var) {
        effect.name <- FF.labs[index.effect, index.pred]
        effect.vals <- FF[, effect.name == pred.names]
        if (any(is.na(effect.vals))) {
          break
        }
        FF[index.effect, index.pred] <- sum(effect.vals * at.mod)
        FF[, index.pred] <- FF[, index.pred] + effect.vals * at.mod[index.effect]

        charge[index.pred, 1] <- charge[index.pred, 1] - at.mod[index.effect] * sum(effect.vals * at.mod)
      }
    }
    new.count.na <- sum(is.na(FF))
    if (count.na == new.count.na) {
      stop("Error: A circularity was detected in the specification of FF. Revise the definition of the linear predictors.")
    }
    count.na <- new.count.na
  }

  list("FF" = FF, "FF.diff" = charge)
}

calc_lin_pred <- function(at, Rt, FF, FF.labs,pred.names) {
  FF.vals <- calc_current_F(at, Rt, FF, FF.labs,pred.names)
  FF <- FF.vals$FF
  FF.diff <- FF.vals$FF.diff
  # print('#########################')
  # print(FF)
  # print(at)
  ft <- crossprod(FF, at) + FF.diff
  Qt <- crossprod(FF, Rt) %*% FF
  list("ft" = ft, "Qt" = Qt, "FF" = FF, "FF.diff" = FF.diff)
}

format_ft <- function(ft, Qt, parms) {
  return(do.call(c, list(ft, Qt)))
}

format_param <- function(conj.param, parms) {
  if (is.null(dim(conj.param))) {
    r.star <- length(conj.param)
    r <- (-1 + sqrt(1 + 4 * r.star)) / 2
    t <- 1
    ft <- conj.param[seq_len(r)] |> matrix(r, t)
    Qt <- conj.param[(r + 1):(r * r + r)] |> array(c(r, r, t))
  } else {
    r.star <- dim(conj.param)[2]
    t <- dim(conj.param)[1]
    r <- (-1 + sqrt(1 + 4 * r.star)) / 2
    ft <- conj.param[, seq_len(r)] |>
      data.frame.to_matrix() |>
      t()
    Qt <- conj.param[, (r + 1):(r * r + r)] |>
      data.frame.to_matrix() |>
      t() |>
      array(c(r, r, t))
  }
  if (t == 1) {
    Qt <- matrix(Qt, r, r)
  }
  return(list("ft" = ft, "Qt" = Qt))
}

generic_param_names <- function(k) {
  index=seq_len(k)
  c(
    paste0("ft.", index),
    paste0(
      "Qt.",
      c(matrix(index, k, k)),
      c(matrix(index, k, k, byrow = TRUE))
    )
  )
}
