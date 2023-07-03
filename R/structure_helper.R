#' Structural blocks for polynomial trends and regressions
#'
#' Creates the structure for a polynomial block with desired order.
#'
#' @param ... Named values for the planning matrix.
#' @param order Positive integer: The order of the polynomial structure.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or scalar: The values for the discount factors associated with the latent variables at each time. If D is an array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and the same discount matrix will be used in all observations. If D is a vector, it should have size t and it is interpreted as the discount factor at each observed time (same discount for all variable). If D is a scalar, the same discount will be used for all latent variables at all times.
#' @param h Matrix, vector or scalar: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x t, where n is the number of latent variables (i.e., the order) and t is the length of the series. If a vector, it should have size T, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.
#' @param H Array, Matrix, vector or scalar: The values for the covariance matrix for the noise factor at each time. If H is an array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the series. If H is a matrix, it's dimensions should be n x n and its values will be used for each time. If H is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of H in the diagonal.
#' @param a1 Vector or scalar: The prior mean for the latent variables associated with this block at time 1. If a1 is a vector, it's dimension should be equal to the order of the polynomial block. If a1 is a scalar, it's value will be used for all latent variables.
#' @param R1 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block at time 1. If R1 is a matrix, it's dimensions should be n x n. If R1 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of R1 in the diagonal.
#' @param monitoring Vector: A vector of flags indicating which variables should be monitored (if automated monitoring is used). Its size should be n. The default is that only the first order component of this structure should be monitored.
#'
#' @return A dlm_block object containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. Its dimension should be n x k x t, where n is the number of latent variables, k is the number of linear predictors in the model and t is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length.
#'    \item h Matrix: The mean for the random noise of the temporal evolution. Its dimension should be n x t.
#'    \item H Array: A 3D-array containing the covariance matrix of the noise for each time. Its dimension should be the same as D.
#'    \item a1 Vector: The prior mean for the latent vector.
#'    \item R1 Matrix: The prior covariance matrix for the latent vector.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (same value as order).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block.
#'    \item monitoring Vector: Same as argument.
#'    \item type Character: The type of block (polynomial).
#' }
#'
#' @export
#' @examples
#' # Creating a first order structure for a model with 2 outcomes.
#' # One block is created for each outcome
#' # with each block being associated with only one of the outcomes.
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 1)
#'
#' # Creating a block with shared effect between the outcomes
#' level_3 <- polynomial_block(alpha1 = 1, alpha2 = 1, order = 2)
#'
#' @details
#'
#' For the ..., D, H, a1 and R1 arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about polynomial trend in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 7.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{search_model}}
#' @family {auxiliary functions for structural blocks}
#'
#'
#' @references
#'    \insertAllCited{}
polynomial_block <- function(..., order = 1, name = "Var_Poly", D = 1, h = 0, H = 0, a1 = 0, R1 = 9, monitoring = c(TRUE, rep(FALSE, order - 1))) {
  if (any(D > 1 | D < 0) & is.numeric(D)) {
    stop("Error: The discount factor D must be a value between 0 and 1 (included).")
  }
  if (any(D == 0)) {
    warning("Some value of D are equal to 0. Those values will be treated as 1.")
  }
  values <- list(...)
  vars <- names(values)
  if (any(vars == "")) {
    stop("Error: One or more linear predictors are unnamed. Please, name all arguments.")
  }
  if (any(vars == "const")) {
    stop("Error: Cannot create a linear predictor named 'const' (reserved name). Choose another label.")
  }
  var_len <- sapply(values, length)
  k <- length(values)
  t <- max(var_len)
  if (any(var_len != t & var_len > 1)) {
    stop(paste0("Error: Outcomes have mismatching lengths. Expected 1 or ", t, ", got: ", paste(var_len, collapse = ", "), "."))
  }

  if (length(D) == 1) {
    D <- array(D, c(order, order, t))
  } else if (is.vector(D)) {
    D <- array(D, c(length(D), order, order)) %>% aperm(c(3, 2, 1))
  } else if (is.matrix(D)) {
    D <- array(D, c(dim(D)[1], dim(D)[2], t))
  }
  t <- if (t == 1) {
    dim(D)[3]
  } else {
    t
  }

  if (length(dim(D)) > 3 | any(dim(D)[1:2] != order) | (dim(D)[3] != t & t > 1)) {
    stop(paste0("Error: Invalid shape for D. Expected ", order, "x", order, "x", t, ". Got ", paste(dim(D), collapse = "x"), "."))
  }

  if (length(H) == 1) {
    pre_H <- diag(order)
    diag(pre_H) <- H
    H <- array(pre_H, c(order, order, t))
  } else if (is.vector(H)) {
    H_vals <- H
    pre_H <- diag(order)
    H <- array(0, c(order, order, length(H_vals)))
    for (i in 1:length(H_vals)) {
      diag(pre_H) <- H_vals[i]
      H[, , i] <- pre_H
    }
  } else if (is.matrix(H)) {
    H <- array(H, c(dim(H)[1], dim(H)[2], t))
  }
  t <- if (t == 1) {
    dim(H)[3]
  } else {
    t
  }

  if (length(dim(h)) < 2) {
    if (t == 1 & length(h) > 1) {
      t <- length(h)
    }
    placeholder <- h
    h <- matrix(0, order, t)
    h[1, ] <- placeholder
  }

  if (any(dim(h) != c(order, t))) {
    stop(paste0("Error: Invalid shape for h. Expected ", order, "x", t, ". Got ", paste(dim(h), collapse = "x"), "."))
  }

  D <- array(D, c(order, order, t))
  H <- array(H, c(order, order, t))

  if (length(dim(H)) > 3 | any(dim(H)[1:2] != order) | (dim(H)[3] != t & t > 1)) {
    stop(paste0("Error: Invalid shape for H. Expected ", order, "x", order, "x", t, ". Got ", paste(dim(H), collapse = "x"), "."))
  }

  FF <- array(0, c(order, k, t), dimnames = list(NULL, vars, NULL))
  FF_labs <- matrix("const", order, k, dimnames = list(NULL, vars))
  for (i in 1:k) {
    name_var <- vars[i]
    if (typeof(values[[name_var]]) == "character") {
      FF[1, i, ] <- NA
      FF_labs[1, i] <- values[[name_var]]
      if (values[[name_var]] == "const") {
        stop("Error: Predictor value is equal to 'const', but 'const' is a reserved name. Choose another label.")
      }
    } else {
      FF[1, i, ] <- values[[name_var]]
    }
  }

  not_observed_flag <- is.na(FF) & (array(FF_labs, c(order, k, t)) == "const")
  D[, , apply(not_observed_flag, 3, any)] <- 1
  H[, , apply(not_observed_flag, 3, any)] <- 0
  FF <- ifelse(not_observed_flag, 0, FF)

  G <- diag(order)
  if (order == 2) {
    G[1, 2] <- 1
  } else if (order > 2) {
    diag(G[1:(order - 1), 2:order]) <- 1
  }

  a1 <- if (length(a1) == 1) {
    rep(a1, order)
  } else {
    a1
  }
  if (length(R1) == 1 | is.vector(R1)) {
    pre_R1 <- diag(order)
    diag(pre_R1) <- R1
    R1 <- pre_R1
  } else {
    R1
  }


  if (length(dim(R1)) > 2) {
    stop(paste0("Error: R1 must be a matrix, but it has ", length(dim(R1)), " dimensions."))
  }
  if (any(dim(R1) != order)) {
    stop(paste0("Error: R1 must have dimensions ", order, "x", order, ". Got ", dim(R1)[1], "x", dim(R1)[2], "."))
  }

  var_labs <- 1:order
  names(var_labs)[1] <- "Level"
  if (order > 1) {
    names(var_labs)[2] <- "Slope"
  }
  if (order >= 3) {
    names(var_labs)[3] <- "Curvature"
  }
  if (order > 3) {
    char_len <- floor(log10(order)) + 1
    names(var_labs)[4:order] <- paste0("Ord_", formatC(3:(order - 1), width = char_len, flag = "0"))
  }

  if (length(monitoring) != order) {
    stop(paste0("Error: monitoring size should be equal to the number of latent variables. Expected ", order, ", got ", length(monitoring), "."))
  }

  var_names <- list()
  var_names[[name]] <- var_labs
  block <- list(
    "FF" = FF,
    "FF_labs" = FF_labs,
    "G" = array(G, c(order, order, t)),
    "G_labs" = matrix("const", order, order),
    "D" = D,
    "h" = h,
    "H" = H,
    "a1" = a1,
    "R1" = R1,
    "var_names" = var_names,
    "order" = order,
    "n" = order,
    "t" = t,
    "k" = k,
    "pred_names" = vars,
    "monitoring" = monitoring,
    "type" = "Polynomial"
  )
  class(block) <- "dlm_block"
  block$status <- check.block.status(block)
  return(block)
}


#' Structural blocks for seasonal trends and regressions
#'
#' Creates the structure for a harmonic block with desired periodicity.
#'
#'
#' @param ... Named values for the planning matrix.
#' @param period Positive integer: The size of the harmonic cycle.
#' @param order Positive integer: The order of the harmonic structure.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or scalar: The values for the discount factors associated with the latent variables at each time. If D is an array, it's dimensions should be (2n) x (2n) x t, where n is the order of the harmonic block and t is the length of the outcomes. If D is a matrix, it's dimensions should be (2n) x (2n) and the same discount matrix will be used in all observations. If D is a vector, it should have size t and it is interpreted as the discount factor at each observed time (same discount for all variable). If D is a scalar, the same discount will be used for all latent variables at all times.
#' @param h Matrix, vector or scalar: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be (2n) x t, where n is the order of the harmonic_block and t is the length of the series. If a vector, it should have size t, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.
#' @param H Array, Matrix, vector or scalar: The values for the covariance matrix for the noise factor at each time. If H is an array, it's dimensions should be (2n) x (2n) x t, where n is the order of the harmonic block and t is the length of the series. If H is a matrix, it's dimensions should be (2n) x (2n) and its values will be used for each time. If H is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of H in the diagonal.
#' @param a1 Vector or scalar: The prior mean for the latent variables associated with this block at time 1. If a1 is a vector, it's dimension should be equal to two times the order of the harmonic block. If a1 is a scalar, it's value will be used for all latent variables.
#' @param R1 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block at time 1. If R1 is a matrix, it's dimensions should be (2n) x (2n). If R1 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of R1 in the diagonal.
#' @param monitoring Vector: A vector of flags indicating which variables should be monitored (if automated monitoring is used). Its size should be 2n. The default is that only the first order component of this structure should be monitored.
#'
#' @return A dlm_block object containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. Its dimension should be n x k x t, where n is the number of latent variables, k is the number of linear predictors in the model and t is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length.
#'    \item h Matrix: The mean for the random noise of the temporal evolution. Its dimension should be n x t.
#'    \item H Array: A 3D-array containing the covariance matrix of the noise for each time. Its dimension should be the same as D.
#'    \item a1 Vector: The prior mean for the latent vector.
#'    \item R1 Matrix: The prior covariance matrix for the latent vector.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item period Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block.
#'    \item monitoring Vector: Same as argument.
#'    \item type Character: The type of block (Harmonic).
#' }
#'
#' @export
#' @examples
#' # Creating seasonal structure for a model with 2 outcomes.
#' # One block is created for each outcome
#' # with each block being associated with only one of the outcomes.
#' season_1 <- harmonic_block(alpha1 = 1, period = 3)
#' season_2 <- harmonic_block(alpha2 = 1, period = 6)
#'
#' # Creating a block with shared effect between the outcomes
#' season_3 <- harmonic_block(alpha = 1, alpha2 = 1, period = 12)
#'
#' @details
#'
#' For the ..., D, H, a1 and R1 arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about the modelling of seasonal trends using harmonics in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 8.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{search_model}}
#' @family {auxiliary functions for structural blocks}
#'
#'
#' @references
#'    \insertAllCited{}
harmonic_block <- function(..., period, order = 1, name = "Var_Sazo", D = 1, h = 0, H = 0, a1 = 0, R1 = 4, monitoring = c(FALSE, FALSE)) {
  w <- 2 * pi / period
  block <- polynomial_block(..., order = 2 * order, name = name, D = D, h = h, H = H, a1 = a1, R1 = R1, monitoring = monitoring)
  G <- matrix(0, 2 * order, 2 * order)
  char_len <- floor(log10(order)) + 1

  G[1:2, 1:2] <- matrix(c(cos(w), -sin(w), sin(w), cos(w)), 2, 2)
  if (order == 1) {
    names(block$var_names[[name]]) <- c("Main", "Aux")
  } else {
    names(block$var_names[[name]])[1:2] <- c(paste0("Main_", formatC(1, width = char_len, flag = "0")), paste0("Aux_", formatC(1, width = char_len, flag = "0")))
    for (i in 2:order) {
      index <- 2 * i - 1
      names(block$var_names[[name]])[index:(index + 1)] <- c(
        paste0("Main_", formatC(i, width = char_len, flag = "0")),
        paste0("Aux_", formatC(i, width = char_len, flag = "0"))
      )

      block$FF[index, , ] <- block$FF[1, , ]
      block$FF_labs[index, ] <- block$FF_labs[1, ]
      G[index:(index + 1), index:(index + 1)] <- matrix(c(cos(w * i), -sin(w * i), sin(w * i), cos(w * i)), 2, 2)
    }
  }

  block$G <- array(G, c(2 * order, 2 * order, block$t))
  block$order <- order
  block$period <- period
  block$type <- "Harmonic"
  return(block)
}

#' Structural blocks for regressions
#'
#' Creates a block for a (dynamic) regression for a covariate X_t.
#'
#' @param ... Named values for the planning matrix.
#' @param lag Non-negative integer: An optional argument providing the maximum lag for the explanatory variables. If a positive value is provided, this block will create additional latent variables to measure the lagged effect of X_t up until the given value. See \insertCite{WestHarr-DLM;textual}{kDGLM}, subsection 9.2.2 item (3).
#' @param zero_fill Bool: A boolean indicating if the block should fill the initial delay values with 0's. If TRUE and lag is positive, the block assumes that X_t=0 for all t<1. If FALSE, the block assumes the useer will provide X_t for all t, such that X_t will have size t+propagation_size
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and its values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param h Matrix, vector or scalar: A drift to be add after the temporal evolution (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be 2 x t, where t is the length of the serie. If a vector, it should have size T, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.
#' @param H Array, Matrix, vector or scalar: The values for the covariance matrix for the noise factor at each time. If H is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If H is a matrix, it's dimensions should be n x n and its values will be used for each time. If H is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of H in the diagonal.
#' @param a1 Vector or scalar: The prior mean for the latent variables associated with this block at time 1. If a1 is a vector, it's dimension should be equal to the order of the polynomial block. If a1 is a scalar, it's value will be used for all latent variables.
#' @param R1 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block at time 1. If R1 is a matrix, it's dimensions should be n x n. If R1 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of R1 in the diagonal.
#' @param monitoring Vector: A vector of flags indicating which variables should be monitored (if automated monitoring is used). Its size should be n. The default is that no variable should be monitored.
#'
#' @return A dlm_block object containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. Its dimension should be n x k x t, where n is the number of latent variables, m is the number of outcomes in the model and t is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. Its dimension should be n x n x t, where n is the number of latent variables and t is the time series length.
#'    \item h Matrix: The mean for the random noise of the temporal evolution. Its dimension should be n x t.
#'    \item H Array: A 3D-array containing the covariance matrix of the noise for each time. Its dimension should be the same as D.
#'    \item a1 Vector: The prior mean for the latent vector.
#'    \item R1 Matrix: The prior covariance matrix for the latent vector.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item lag Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block.
#'    \item monitoring Vector: Same as argument.
#'    \item type Character: The type of block (Harmonic).
#' }
#'
#' @export
#' @examples
#'
#' T <- 200
#' X <- rgamma(T, 5, 5)
#' data <- rpois(T, exp(2 + 1.5 * X))
#'
#' level <- polynomial_block(rate = 1, D = 0.95)
#' regression <- regression_block(rate = X, D = 0.98)
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, regression, outcomes = outcome)
#' summary(fitted_data)
#'
#' @details
#'
#' For the ..., D, H, a1 and R1 arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{search_model}}
#' @family {auxiliary functions for structural blocks}
#'
#'
#' @references
#'    \insertAllCited{}
regression_block <- function(..., lag = 0, zero_fill = TRUE, name = "Var_Reg", D = 1, h = 0, H = 0, a1 = 0, R1 = 9, monitoring = rep(FALSE, lag + 1)) {
  order <- lag + 1
  block <- polynomial_block(..., order = order, name = name, D = D, h = h, H = H, a1 = a1, R1 = R1, monitoring = monitoring)

  G <- diag(lag + 1)
  block$G <- array(G, c(order, order, block$t))
  block$order <- NULL
  block$lag <- lag
  block$zero_fill <- zero_fill
  block$type <- "Regression"

  names(block$var_names[[name]]) <- paste0("Lag_", 0:(order - 1))

  t <- block$t
  if (lag > 0) {
    for (i in 1:lag) {
      block$FF[i + 1, , -1] <- block$FF[i, , -t]
    }
  }
  if (!zero_fill & t > 1) {
    block$FF <- block$FF[, , (lag + 1):t, drop = FALSE]
    block$G <- block$G[, , (lag + 1):t, drop = FALSE]
    block$D <- block$D[, , (lag + 1):t, drop = FALSE]
    block$h <- block$h[, (lag + 1):t, drop = FALSE]
    block$H <- block$H[, , (lag + 1):t, drop = FALSE]
    block$t <- t - lag
  }

  return(block)
}

#' Structural blocks for auto regressive trends and regressions
#'
#' Creates the structure for a Auto Regressive (AR) block (see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 9) with desired order.
#' As the package suppose that the structure of the model is linear, a linearization is applied to the evolution equation, as described in \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 13.
#' This block also supports Transfer Functions, being necessary to specify the associated pulse when calling the AR_block function (see arg.).
#'
#' @param ... Named values for the planning matrix.
#' @param order Positive integer: The order of the AR block.
#' @param noise_var Non-negative scalar: The variance of the white noise added to the latent state.
#' @param noise_discount Vector or scalar: The value for the discount factor associated with the current latent state. If noise_discount is a vector, it should have size t and it is interpreted as the discount factor at each observed time. If D is a scalar, the same discount will be used for all observation.
#' @param pulse Vector or scalar: An optional argument providing the values for the pulse for a Transfer Function. Default is 0 (no Transfer Function).
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param AR_support String: Either "constrained" or "free" (default). If AR_support is "constrained", then the AR coefficients will be forced to be on the interval (-1,1), otherwise, the coefficients will be unrestricted. Beware that, under no restriction on the coefficients, there is no guarantee that the estimated coefficients will imply in a stationary process, furthermore, if the order of the AR block is greater than 1. As such the restriction of the coefficients support is only available for AR blocks with order equal to 1.
#' @param a1 Vector or scalar: The prior mean for the AR coefficients associated with this block at time 1. If a1 is a vector, it's dimension should be equal to the order of the AR block. If a1 is a scalar, it's value will be used for all coefficients. If the coefficients are restricted to the interval (-1,1), the a1 is interpreted as the mean for atanh(rho), where rho is the AR coefficient.
#' @param R1 Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with this block at time 1. If R1 is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If R1 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of R1 in the diagonal. If the coefficients are restricted to the interval (-1,1), the R1 is interpreted as the covariance matrix for atanh(rho), where rho is the AR coefficient.
#' @param monitoring Vector: A vector of flags indicating which AR coefficients should be monitored (if automated monitoring is used). Its size should be n, where n is the order of the AR block. The default is that no coefficient should be monitored.
#' @param D Array, Matrix, vector or scalar: The values for the discount factors associated with the AR coefficients at each time. If D is an array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and the same discount matrix will be used in all observations. If D is a vector, it should have size t and it is interpreted as the discount factor at each observed time (same discount for all variable). If D is a scalar, the same discount will be used for all AR coefficients at all times.
#' @param h Matrix, vector or scalar: A drift to be add in the AR coefficients after the temporal evolution (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the order of the AR block and T is the length of the series. If a scalar, the passed value will be used for all coefficients at each time.
#' @param H Array, Matrix, vector or scalar: The values for the covariance matrix for the noise factor associated with the AR coefficients at each time. If H is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If H is a matrix, it's dimensions should be n x n and its values will be used for each time. If H is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of H in the diagonal.
#' @param a1_states Vector or scalar: The prior mean for the states associated with this block at time 1. If a1_states is a vector, it's dimension should be equal to the order of the AR block. If a1_states is a scalar, it's value will be used for all coefficients.
#' @param R1_states Matrix, vector or scalar: The prior covariance matrix for the states associated with this block at time 1. If R1_states is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If R1_state is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of R1_state in the diagonal.
#' @param monitoring_states bool: A flag indicating if the latent state should be monitored (if automated monitoring is used). The default is TRUE.
#' @param h_states Vector or scalar: A drift to be add in the states after the temporal evolution (can be interpreted as the mean of the random noise at each time). If a vector, it should have size t, and each value will be applied in their respective time. If a scalar, the passed value will be used for all observations.
#' @param a1_pulse Vector or scalar: The prior mean for the coefficients associated with the pulses at time 1. If a1_pulse is a vector, it's dimension should be equal to the number of pulses. If a1_pulse is a scalar, it's value will be used for all coefficients.
#' @param R1_pulse Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with the pulses at time 1. If R1_pulse is a matrix, it's dimensions should be n x n, where n is the number of pulses. If R1_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of R1_pulse in the diagonal.
#' @param D_pulse Array, Matrix, vector or scalar: The values for the discount factors associated with the pulse coefficients at each time. If D_pulse is an array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If D_pulse is a matrix, it's dimensions should be n x n and the same discount matrix will be used in all observations. If D_pulse is a vector, it should have size t and it is interpreted as the discount factor at each observed time (same discount for all variable). If D is a scalar, the same discount will be used for all pulse coefficients at all times.
#' @param h_pulse Matrix, vector or scalar: A drift to be add in the pulse effect after the temporal evolution (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the number of pulses and T is the length of the series. If a scalar, the passed value will be used for all latent variable at each time.
#' @param H_pulse Array, Matrix, vector or scalar: The values for the covariance matrix for the noise factor associated with pulse coefficients at each time. If H_pulse is an array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If H_pulse is a matrix, it's dimensions should be n x n and its values will be used for each time. If H_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of H_pulse in the diagonal.
#' @param monitoring_pulse Vector: A vector of flags indicating which pulse coefficients should be monitored (if automated monitoring is used). Its size should be n, where n is the number of pulses. The default is that no pulse coefficient should be monitored.
#'
#' @return A dlm_block object containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item H Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item a1 Vector: The prior mean for the latent vector.
#'    \item R1 Matrix: The prior covariance matrix for the latent vector.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block.
#'    \item monitoring Vector: The combination of monitoring, monitoring_states and monitoring_pulse.
#'    \item type Character: The type of block (AR).
#'    \item AR_support Character: Same as argument.
#' }
#'
#' @export
#' @examples
#'
#' AR_block(mu = 1, pulse = rnorm(200), order = 3, noise_var = 0.1)
#'
#' @details
#'
#' For the ..., noise_var, D, H, a1, R1, a1_states, R1_states, D_states, a1_pulse, R1_pulse, D_pulse, S_pulse arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about Auto regressive models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 9.
#'
#' For the details about the linearization of non-linear evolution equations in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 13.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @family {auxiliary functions for structural blocks}
#'
#' @references
#'    \insertAllCited{}
AR_block <- function(..., order, noise_var, noise_discount = 1, pulse = 0, name = "Var_AR", AR_support = "free",
                     D = 1, h = 0, H = 0, a1 = c(0, rep(0, order - 1)), R1 = 0.25, monitoring = rep(FALSE, order),
                     a1_states = 0, R1_states = 9, h_states = 0, monitoring_states = TRUE,
                     a1_pulse = 0, R1_pulse = 9, D_pulse = 1, h_pulse = 0, H_pulse = 0, monitoring_pulse = NA) {
  if (AR_support == "constrained" & order > 1) {
    warning("Restrictions on the AR coefficients are only available for first order models. Ignoring AR_support argument.")
    AR_support <- "free"
  }

  h_states_placeholder <- h_states
  h_states <- matrix(0, order, length(h_states))
  h_states[1, ] <- h_states_placeholder
  H_states <- array(0, c(order, order, length(noise_var)))
  H_states[1, 1, ] <- noise_var
  D_states <- array(1, c(order, order, length(noise_var)))
  D_states[1, 1, ] <- noise_discount
  block_state <-
    polynomial_block(..., order = order, name = paste0(name, "_State"), a1 = a1_states, R1 = R1_states, D = D_states, h = h_states, H = H_states, monitoring = c(monitoring_states, rep(FALSE, order - 1)))


  names(block_state$var_names[[paste0(name, "_State")]]) <- paste0("Lag_", 0:(order - 1))

  dummy_var <- list()
  dummy_var[[names(list(...))[1]]] <- rep(0, block_state$t)
  if (length(h_states) > 1 & length(dim(h_states)) < 2) {
    h_states <- matrix(h_states, order, length(h_states))
  }
  block_coeff <-
    do.call(
      function(...) {
        polynomial_block(..., order = order, name = paste0(name, "_Coef"), a1 = a1, R1 = R1, D = D, h = h_states, H = H, monitoring = monitoring)
      },
      dummy_var
    )


  names(block_coeff$var_names[[paste0(name, "_Coef")]]) <- paste0("Lag_", 0:(order - 1))

  block <- block_state + block_coeff
  k <- block$k

  if (order == 1) {
    G <- diag(2)
    G_labs <- matrix("const", 2, 2)
    G[1, 1] <- NA
    G_labs[1, 1] <- tolower(AR_support)
  } else {
    G <- matrix(0, 2 * order, 2 * order)
    G_labs <- matrix("const", 2 * order, 2 * order)
    G[1, 1:order] <- NA
    G_labs[1, 2 * (1:order) - 1] <- tolower(AR_support)
    G[2:order, -(order:(2 * order))] <- diag(order - 1)
    G[(order + 1):(2 * order), (order + 1):(2 * order)] <- diag(order)
    index <- sort(c(c(1:order), c(1:order)))
    G <- G[index + c(0, order), index + c(0, order)]
  }
  block$G <- array(G, c(2 * order, 2 * order, block$t))
  block$G_labs <- G_labs
  block$order <- order
  block$type <- "AR"
  block$AR_support <- AR_support
  true_order <- c(rbind(1:order, 1:order + order))
  block$a1 <- block$a1[true_order]
  block$R1 <- block$R1[true_order, true_order]
  block$D <- block$D[true_order, true_order, ]
  block$H <- block$H[true_order, true_order, ]
  if (any(pulse != 0)) {
    k <- if.null(dim(pulse)[2], 1)
    t <- if.null(dim(pulse)[1], length(pulse))
    dummy_var <- list()
    dummy_var[[names(list(...))[1]]] <- rep(0, t)
    if (length(h_pulse) > 1 & length(dim(h_pulse)) < 2) {
      h_pulse <- matrix(h_pulse, order, length(h_pulse))
    }
    if (is.na(monitoring_pulse)) {
      monitoring_pulse <- rep(FALSE, k)
    }
    block_pulse <-
      do.call(
        function(...) {
          polynomial_block(..., order = k, name = paste0(name, "_Pulse"), a1 = a1_pulse, R1 = R1_pulse, D = D_pulse, h = h_pulse, H = H_pulse, monitoring = monitoring_pulse)
        },
        dummy_var
      )

    names(block_pulse$var_names[[paste0(name, "_Pulse")]]) <- paste0("Pulse_", 1:k)

    block_pulse$G <- diag(k)
    block <- block + block_pulse
    block$G[1, (2 * order + 1):(2 * order + k), ] <- t(matrix(pulse, block$t, k))
  }
  return(block)
}

#' Auxiliary function for block superposition
#'
#' An auxiliary function for block superposition.
#'
#' @param ... dlm_block: A sequence of block to be combine.
#'
#' @return The combined blocks as a dlm_block.
#' @export
#'
#' @examples
#'
#' # Long way
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 2)
#' season_2 <- harmonic_block(alpha2 = 1, period = 20)
#'
#' final_block <- block_superpos(level_1, level_2, season_2)
#'
#' # Short way
#' final_block <- polynomial_block(alpha1 = 1, order = 1) +
#'   polynomial_block(alpha2 = 1, order = 2) +
#'   harmonic_block(alpha2 = 1, period = 20)
#'
#' @details
#' Additional details can be found in \insertCite{WestHarr-DLM;textual}{kDGLM}, section 6.2.
#'
#'
#' @family {auxiliary functions for structural blocks}
#'
#' @references
#'    \insertAllCited{}
block_superpos <- function(...) {
  blocks <- list(...)
  if (length(blocks) == 1) {
    return(blocks[[1]])
  }
  for (block in blocks) {
    if (!inherits(block, "dlm_block")) {
      stop(paste0("Error: Expected all arguments to be dlm_block's, but got a ", class(block), "."))
    }
  }

  n <- 0
  t <- 1
  k <- 1

  var_names <- list()
  names <- c()
  pred_names <- c()
  count_labs <- 0
  for (block in blocks) {
    ref_names <- block$var_names
    count_labs_i <- 0
    for (name in names(ref_names)) {
      count_labs <- count_labs + 1
      count_labs_i <- count_labs_i + 1
      names <- c(names, name)
      var_names[[count_labs]] <- ref_names[[count_labs_i]] + n
    }
    pred_names <- c(pred_names, block$pred_names)
    if (block$t > 1) {
      if (block$t != t & t > 1) {
        stop(paste("Error: Blocks should have same length or length equal 1. Got", block$t, "and", t))
      }
      t <- block$t
    }
    n <- n + block$n
  }
  names(var_names) <- names

  # for(name in unique(names)){
  #   count_name=sum(names==name)
  #   if(count_name>1){
  #     len_char=floor(log10(count_name))+1
  #     names[names==name]=paste0(name,'_',formatC(1:count_name,width=len_char,flag='0'))
  #   }
  # }

  pred_names <- sort(unique(pred_names))
  k <- length(pred_names)

  FF <- array(0, c(n, k, t), dimnames = list(NULL, pred_names, NULL))
  FF_labs <- matrix("const", n, k, dimnames = list(NULL, pred_names))
  G <- array(0, c(n, n, t))
  G_labs <- matrix("const", n, n)
  D <- array(0, c(n, n, t))
  h <- matrix(0, n, t)
  H <- array(0, c(n, n, t))
  a1 <- c()
  R1 <- matrix(0, n, n)
  position <- 1
  monitoring <- c()
  status <- "defined"
  for (block in blocks) {
    k_i <- length(block$pred_names)
    current_range <- position:(position + block$n - 1)
    FF[current_range, pred_names %in% block$pred_names, ] <- block$FF[, (1:k_i)[order(block$pred_names)], ]
    FF_labs[current_range, pred_names %in% block$pred_names] <- block$FF_labs[, (1:k_i)[order(block$pred_names)]]

    G[current_range, current_range, ] <- block$G
    G_labs[current_range, current_range] <- block$G_labs
    D[current_range, current_range, ] <- block$D
    h[current_range, ] <- block$h
    H[current_range, current_range, ] <- block$H
    a1 <- c(a1, block$a1)
    R1[current_range, current_range] <- block$R1
    position <- position + block$n
    monitoring <- c(monitoring, block$monitoring)
  }
  block <- list(
    "FF" = FF,
    "FF_labs" = FF_labs,
    "G" = G,
    "G_labs" = G_labs,
    "D" = D,
    "h" = h,
    "H" = H,
    "a1" = a1,
    "R1" = R1,
    "n" = n,
    "t" = t,
    "k" = k,
    "status" = status,
    "var_names" = var_names,
    "name" = names,
    "pred_names" = pred_names,
    "monitoring" = monitoring,
    "type" = "Mixed"
  )
  class(block) <- "dlm_block"
  block$status <- check.block.status(block)
  return(block)
}

#' Auxiliary function to replicate blocks
#'
#' An auxiliary function to replicate blocks.
#'
#' @param block dlm_block: A block to be replicated
#' @param k Integer: The number of blocks to generate.
#'
#' @return The combined replicated blocks as a dlm_block.
#' @export
#'
#' @examples
#' # Long way
#' level <- polynomial_block(alpha = 1, order = 1)
#'
#' final_block <- block_mult(level, 5)
#'
#' # Short way
#' final_block <- 5 * polynomial_block(alpha = 1, order = 1)
#'
#' @seealso \code{\link{block_rename}}
#' @family {auxiliary functions for structural blocks}
block_mult <- function(block, k) {
  block_list <- list()
  size_total <- floor(log10(k)) + 1
  if (k > 1) {
    block_ref <- block
    block$pred_names <- paste0(block$pred_names, "_", paste0(rep("0", size_total - 1), collapse = ""), "1")
    block_list[[1]] <- block
    for (i in 2:k) {
      size_i <- floor(log10(i)) + 1
      block_clone <- block_ref
      block_clone$pred_names <- paste0(block_ref$pred_names, "_", paste0(rep("0", size_total - size_i), collapse = ""), i)
      block_list[[i]] <- block_clone
    }
    block <- do.call(block_superpos, block_list)
  }
  return(block)
}

#' block_rename
#'
#' @param block A dlm_block object.
#' @param pred_names A vector of string with names for each linear predictor in block.
#'
#' @return A dlm_block with the linear predictors renamed to the values passed in names.
#' @export
#'
#' @examples
#'
#' base_block <- polynomial_block(
#'   eta = 1,
#'   order = 1,
#'   name = "Poly",
#'   D = 0.95
#' )
#'
#' final_block <- block_rename(2 * base_block, c("mu", "sigma"))
#'
#' @family {auxiliary functions for structural blocks}
block_rename <- function(block, pred_names) {
  if (!inherits(block, "dlm_block")) {
    stop("Error: The block argument is not a dlm_block object.")
  }
  if (length(pred_names) != length(block$pred_names)) {
    stop(paste0("Error: The number of names provided does not match the number of linear predictor in the block. Expected ", length(block$pred_names), ", got ", length(names), "."))
  }
  if (length(pred_names) != length(unique(pred_names))) {
    stop(paste0("Error: Repeated names are not allowed."))
  }

  block$pred_names <- pred_names
  colnames(block$FF) <- pred_names
  block$status <- check.block.status(block)
  return(block)
}
