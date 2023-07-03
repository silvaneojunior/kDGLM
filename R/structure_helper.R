#' Structural blocks for polynomial trends and regressions
#'
#' Creates the structure for a polynomial block with desired order.
#'
#' @param ... Named values for the planing matrix.
#' @param order Positive integer: The order of the polynomial structure.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimension should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimensions should be n x n. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the first element of C0 is equal to NA, then a proper variance is calculated for the first latent state based on the scale of the effect of that state on the linear predictor.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containing the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (same value as order).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
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
#' For the ..., D, W, m0 and C0 arguments, the user may set one or more of it's values as a string.
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
polynomial_block <- function(..., order = 1, name = "Var_Poly", D = 1, W = 0, m0 = 0, C0 = c(NA, rep(1, order - 1))) {
  if (any(D > 1 | D < 0) & is.numeric(D)) {
    stop("Error: The discount factor D must be a value between 0 and 1 (included).")
  }
  if (any(D == 0)) {
    warning("Some value of D are equal to 0. Those values will be treated as 1.")
  }
  values <- list(...)
  vars <- names(values)
  if (any(vars == "")) {
    stop("Error: One or more outcomes are unnamed. Please, name all arguments.")
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

  if (length(W) == 1) {
    pre_W <- diag(order)
    diag(pre_W) <- W
    W <- array(pre_W, c(order, order, t))
  } else if (is.vector(W)) {
    W_vals <- W
    pre_W <- diag(order)
    W <- array(0, c(order, order, length(W_vals)))
    for (i in 1:length(W_vals)) {
      diag(pre_W) <- W_vals[i]
      W[, , i] <- pre_W
    }
  } else if (is.matrix(W)) {
    W <- array(W, c(dim(W)[1], dim(W)[2], t))
  }
  t <- if (t == 1) {
    dim(W)[3]
  } else {
    t
  }
  D <- array(D, c(order, order, t))

  if (length(dim(W)) > 3 | any(dim(W)[1:2] != order) | (dim(W)[3] != t & t > 1)) {
    stop(paste0("Error: Invalid shape for W. Expected ", order, "x", order, "x", t, ". Got ", paste(dim(W), collapse = "x"), "."))
  }

  for (name_var in vars[var_len == 1]) {
    var_len[[name_var]] <- t
    values[[name_var]] <- rep(values[[name_var]], t)
  }

  FF <- array(0, c(order, k, t))
  FF[1, , ] <- matrix(sapply(values, c), k, t, byrow = TRUE)
  D[, , apply(is.na(FF), 3, any)] <- 1
  W[, , apply(is.na(FF), 3, any)] <- 0
  FF <- ifelse(is.na(FF), 0, FF)

  G <- diag(order)
  if (order == 2) {
    G[1, 2] <- 1
  } else if (order > 2) {
    diag(G[1:(order - 1), 2:order]) <- 1
  }

  m0 <- if (length(m0) == 1) {
    rep(m0, order)
  } else {
    m0
  }
  if (length(C0) == 1 | is.vector(C0)) {
    pre_C0 <- diag(order)
    diag(pre_C0) <- C0
    C0 <- pre_C0
  } else {
    C0
  }
  if (any(is.na(C0))) {
    ref_val <- as.numeric(c(...))
    ref_val <- ref_val[!is.na(ref_val)]
    if (all(ref_val == 0 | ref_val == 1) | if.na(var(ref_val), 0) == 0) {
      ref_var <- 1
    } else {
      ref_var <- 1 / var(ref_val)
    }
    if (is.na(C0[1, 1])) {
      C0[1, 1] <- ref_var
    }
    C0[is.na(C0)] <- 1
  }


  if (length(dim(C0)) > 2) {
    stop(paste0("Error: C0 must be a matrix, but it has ", length(dim(C0)), " dimensions."))
  }
  if (any(dim(C0) != order)) {
    stop(paste0("Error: C0 must have dimensions ", order, "x", order, ". Got ", dim(C0)[1], "x", dim(C0)[2], "."))
  }

  names <- list()
  names[[name]] <- c(1:order)
  block <- list(
    "FF" = FF,
    "G" = array(G, c(order, order, t)),
    "G_labs" = matrix("const", order, order),
    "D" = D,
    "W" = W,
    "m0" = m0,
    "C0" = C0,
    "names" = names,
    "order" = order,
    "n" = order,
    "t" = t,
    "k" = k,
    "var_names" = vars,
    "type" = "Polynomial"
  )
  class(block) <- "dlm_block"
  if (any(is.character(if.na(FF, 0))) |
    any(is.character(if.na(G, 0))) |
    any(is.character(if.na(D, 0))) |
    any(is.character(if.na(W, 0))) |
    any(is.character(if.na(m0, 0))) |
    any(is.character(if.na(C0, 0)))
  ) {
    block$status <- "undefined"
  } else {
    block$status <- "defined"
  }
  return(block)
}


#' Structural blocks for seasonal trends and regressions
#'
#' Creates the structure for a harmonic block with desired periodicity.
#'
#' @param ... Named values for the planing matrix.
#' @param period Positive integer: The size of the harmonic cycle.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimension should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimensions should be n x n. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containing the variables indexes by their name.
#'    \item period Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
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
#' For the ..., D, W, m0 and C0 arguments, the user may set one or more of it's values as a string.
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
harmonic_block <- function(..., period, name = "Var_Sazo", D = 1, W = 0, m0 = 0, C0=c(NA, 1)) {
  w <- 2 * pi / period
  order <- 2
  block <- polynomial_block(..., order = order, name = name, D = D, W = W, m0 = m0, C0 = C0)

  G <- matrix(c(cos(w), -sin(w), sin(w), cos(w)), order, order)
  block$G <- array(G, c(order, order, block$t))
  block$order <- NULL
  block$period <- period
  block$type <- "Harmonic"
  return(block)
}

#' Structural blocks for auto regressive trends and regressions
#'
#' Creates the structure for a Auto Regressive (AR) block (see West & Harrison (1997), chapter 9) with desired order.
#' As the package suppose that the structure of the model is linear, a linearization is applied to the evolution equation, as described in West & Harrison (1997), chapter 13.
#' This block also supports Transfer Functions, being necessary to specify the associated pulse when calling the AR_block function (see arg.).
#'
#' @param ... Named values for the planing matrix.
#' @param order Positive integer: The order of the AR block.
#' @param noise_var Non negative scalar: The variance of the white noise added to the latent state.
#' @param noise_var Non negative scalar: The variance of the white noise added to the latent state.
#' @param pulse Vector or scalar: An optional argument providing the values for the pulse for a Transfer Function. Default is 0.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param AR_support String: Either "constrained" or "free". If AR_support is "constrained", then the AR coefficients will be forced to be on the interval (-1,1), otherwise, the coefficients will be unrestricted. Beware that, under no restriction on the coefficients, there is no guarantee that the estimated coefficients will imply in a stationary process, furthermore, if the order of the AR block is greater than 1, then the restriction imposed when AR_support is equal to "constrained" does NOT guarantee that the process will be stationary (although it may help).
#' @param m0 Vector or scalar: The prior mean for the coefficients associated with this block. If m0 is a vector, it's dimension should be equal to the order of the AR block. If m0 is a scalar, it's value will be used for all coefficients. If the coefficients are restricted to the interval (-1,1), the m0 is interpreted as the mean for logit((rho+1)/2), where rho is the AR coefficient.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with this block. If C0 is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the coefficients are restricted to the interval (-1,1), the C0 is interpreted as the covariance matrix for logit((rho+1)/2), where rho is the AR coefficient.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors associated with the AR coefficients at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor associated with the AR coefficients at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the coefficients associated with this block. If m0 is a vector, it's dimension should be equal to the order of the AR block. If m0 is a scalar, it's value will be used for all coefficients. If the coefficients are restricted to the interval (-1,1), the m0 is interpreted as the mean for logit((rho+1)/2), where rho is the AR coefficient.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with this block. If C0 is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the coefficients are restricted to the interval (-1,1), the C0 is interpreted as the covariance matrix for logit((rho+1)/2), where rho is the AR coefficient.
#' @param m0_states Vector or scalar: The prior mean for the states associated with this block. If m0_states is a vector, it's dimension should be equal to the order of the AR block. If m0_states is a scalar, it's value will be used for all coefficients.
#' @param C0_states Matrix, vector or scalar: The prior covariance matrix for the states associated with this block. If C0_states is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0_state is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0_state in the diagonal.
#' @param D_states Array, Matrix, vector or  scalar: The values for the discount factors for the states associated with this block. If D_states is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D_states is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D_states is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D_states in the diagonal.
#' @param m0_pulse Vector or scalar: The prior mean for the coefficients associated with the pulses. If m0_pulse is a vector, it's dimension should be equal to the number of pulses. If m0_pulse is a scalar, it's value will be used for all coefficients.
#' @param C0_pulse  Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with the pulses. If C0_pulse is a matrix, it's dimensions should be n x n, where n is the number of pulses. If C0_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0_pulse in the diagonal.
#' @param D_pulse Array, Matrix, vector or  scalar: The values for the discount factors associated with pulse coefficients at each time. If D_pulse is a array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If D_pulse is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D_pulse is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D_pulse in the diagonal.
#' @param W_pulse Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor associated with pulse coefficients at each time. If W_pulse is a array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If W_pulse is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of W_pulse in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containing the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
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
#' For the ..., noise_var, D, W, m0, C0, m0_states, C0_states, D_states, m0_pulse, C0_pulse, D_pulse, W_pulse arguments, the user may set one or more of it's values as a string.
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
AR_block <- function(..., order, noise_var, pulse = 0, name = "Var_AR", AR_support = "free",
                     D = 1, W = 0, m0 = c(1,rep(0,order-1)), C0 = 1,
                     m0_states = 0, C0_states = c(NA, rep(1, order - 1)), D_states = 1,
                     m0_pulse = 0, C0_pulse = 1, D_pulse = 1, W_pulse = 0) {
  W_states=diag(order)*0
  W_states[1,1]=noise_var
  block_state <-
    polynomial_block(..., order = order, name = paste0(name, "_State"), m0 = m0_states, C0 = C0_states, D = D_states, W = W_states)


  dummy_var <- list()
  dummy_var[[names(list(...))[1]]] <- rep(0, block_state$t)
  block_coeff <-
    do.call(
      function(...) {
        polynomial_block(..., order = order, name = paste0(name, "_Coeff"), m0 = m0, C0 = C0, D = D, W = W)
      },
      dummy_var
    )

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
    G_labs[1, 2*(1:order)-1] <- tolower(AR_support)
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
  block$m0 <- block$m0[true_order]
  block$C0 <- block$C0[true_order, true_order]
  block$D <- block$D[true_order, true_order, ]
  block$W <- block$W[true_order, true_order, ]
  block$names[[1]] <- seq(1, 2 * order, 2)[block$names[[1]]]
  block$names[[2]] <- seq(2, 2 * order, 2)[block$names[[2]] - order]
  if (any(pulse != 0)) {
    k <- if.null(dim(pulse)[2], 1)
    t <- if.null(dim(pulse)[1], length(pulse))
    dummy_var <- list()
    dummy_var[[names(list(...))[1]]] <- rep(0, t)
    block_pulse <-
      do.call(
        function(...) {
          polynomial_block(..., order = k, name = paste0(name, "_Pulse"), m0 = m0_pulse, C0 = C0_pulse, D = D_pulse, W = W_pulse)
        },
        dummy_var
      )
    block_pulse$G <- diag(k)
    block <- block + block_pulse
    block$G[1, (2 * order + 1):(2 * order + k), ] <- t(matrix(pulse, block$t, k))
  }
  return(block)
}

#' correlation_block
#'
#' DESCRIPTION
#'
#' @param var_names Vector: Name of the linear predictors associated with this block.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimension should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimensions should be n x n. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containing the variables indexes by their name.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (Correlation).
#' }
#'
#' @export
#' @keywords internal
#' @examples
#' # EXAMPLE
#'
#' @family {auxiliary functions for structural blocks}
#'
#'
#' @references
#'    \insertAllCited{}
correlation_block <- function(var_names, order = 1, name = "Var_cor", D = 1, W = 0, m0 = 0, C0 = 1) {
  k <- length(var_names)
  arg_list <- c(rep(0, k), list(order = k, name = name, D = D, W = W, m0 = m0, C0 = C0))
  names(arg_list) <- c(var_names, "order", "name", "D", "W", "m0", "C0")
  block <- do.call(polynomial_block, arg_list)

  block$FF <- simplify2array(apply(block$FF, 3, function(x) {
    diag(x) <- NA
    rbind(x, NA)
  },
  simplify = FALSE
  ))

  W <- 1
  block$G <- array(diag(c(rep(1, k), 0)), c(k + 1, k + 1, block$t))
  block$D <- simplify2array(apply(block$D, 3, function(x) {
    as.matrix(bdiag(x, 1))
  },
  simplify = FALSE
  ))
  block$W <- simplify2array(apply(block$W, 3, function(x) {
    as.matrix(bdiag(x, W))
  },
  simplify = FALSE
  ))
  block$m0 <- c(block$m0, 0)
  block$C0 <- as.matrix(bdiag(block$C0, W))
  block$k <- k

  n <- k + 1
  block$names[[name]] <- c(block$names[[name]], n)
  block$order <- order
  block$n <- n
  if (order > 1) {
    block_ref <- block
    for (i in 2:order) {
      block_copy <- block_ref
      block_copy$m0[1:(i - 1)] <- 0
      block_copy$C0[1:(i - 1), ] <- 0
      block_copy$C0[, 1:(i - 1)] <- 0
      block <- block + block_copy
    }
  }


  block$type <- "Correlation"
  return(block)
}

#' Auxiliary function to merge blocks
#'
#' An auxiliary function to merge blocks.
#'
#' @param ... dlm_block: A sequence of block to be merged.
#'
#' @return The merged block as a dlm_block.
#' @export
#'
#' @examples
#'
#' # Long way
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 2)
#' season_2 <- harmonic_block(alpha2 = 1, period = 20)
#'
#' final_block <- block_merge(level_1, level_2, season_2)
#'
#' # Short way
#' final_block <- polynomial_block(alpha1 = 1, order = 1) +
#'   polynomial_block(alpha2 = 1, order = 2) +
#'   harmonic_block(alpha2 = 1, period = 20)
#'
#' @family {auxiliary functions for structural blocks}
block_merge <- function(...) {
  blocks <- list(...)
  for (block in blocks) {
    if (!inherits(block, "dlm_block")) {
      stop(paste0("Error: Expected all arguments to be dlm_block's, but got a ", class(block), "."))
    }
  }

  n <- 0
  t <- 1
  k <- 1

  names <- list()
  var_names <- c()
  for (block in blocks) {
    ref_names <- block$names
    for (name in names(ref_names)) {
      ref_names[[name]] <- ref_names[[name]] + n
      names[[name]] <- c(names[[name]],ref_names[[name]])
    }
    # names <- c(names, ref_names)
    var_names <- c(var_names, block$var_names)
    if (block$t > 1) {
      if (block$t != t & t > 1) {
        stop(paste("Error: Blocks should have same length or length equal 1. Got", block$t, "and", t))
      }
      t <- block$t
    }
    n <- n + block$n
  }
  var_names <- sort(unique(var_names))
  k <- length(var_names)
  # for (name in names(names)) {
  #   ref_idx <- which(names(names) == name)
  #   n_names <- length(ref_idx)
  #   if (n_names > 1) {
  #     names(names)[ref_idx] <- paste0(names(names)[ref_idx], "_", c(1:n_names))
  #   }
  # }

  FF <- array(0, c(n, k, t), dimnames = list(NULL, var_names, NULL))
  G <- array(0, c(n, n, t))
  G_labs <- matrix("const", n, n)
  D <- array(0, c(n, n, t))
  W <- array(0, c(n, n, t))
  m0 <- c()
  C0 <- matrix(0, n, n)
  position <- 1
  status <- "defined"
  for (block in blocks) {
    k_i <- length(block$var_names)
    current_range <- position:(position + block$n - 1)
    FF[current_range, var_names %in% block$var_names, ] <- block$FF[, (1:k_i)[order(block$var_names)], ]

    G[current_range, current_range, ] <- block$G
    G_labs[current_range, current_range] <- block$G_labs
    D[current_range, current_range, ] <- block$D
    W[current_range, current_range, ] <- block$W
    m0 <- c(m0, block$m0)
    C0[current_range, current_range] <- block$C0
    position <- position + block$n
    status <- if (block$status == "undefined") {
      "undefined"
    } else {
      status
    }
  }
  block <- list(
    "FF" = FF,
    "G" = G,
    "G_labs" = G_labs,
    "D" = D,
    "W" = W,
    "m0" = m0,
    "C0" = C0,
    "n" = n,
    "t" = t,
    "k" = k,
    "status" = status,
    "names" = names,
    "var_names" = var_names
  )
  class(block) <- "dlm_block"
  return(block)
}

#' Auxiliary function to replicate blocks
#'
#' An auxiliary function to merge blocks.
#'
#' @param block dlm_block: A block to be multiplied.
#' @param k Integer: The number of blocks to generate.
#'
#' @return The merged block as a dlm_block.
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
  size_total=floor(log10(k))+1
  if (k > 1) {
    block_ref <- block
    block$var_names <- paste0(block$var_names, "_",paste0(rep('0',size_total-1),collapse=''),"1")
    block_list[[1]] <- block
    for (i in 2:k) {
      size_i=floor(log10(i))+1
      block_clone <- block_ref
      block_clone$var_names <- paste0(block_ref$var_names, "_",paste0(rep('0',size_total-size_i),collapse=''), i)
      block_list[[i]] <- block_clone
    }
    block <- do.call(block_merge, block_list)
  }
  return(block)
}

#' block_rename
#'
#' @param block A dlm_block object.
#' @param names A vector of string with names for each linear predictor in block.
#'
#' @return A dlm_block with the linear predictors renamed to the values passed in names.
#' @export
#'
#' @examples
#'
#' base_block=polynomial_block(
#' eta=1,
#' order = 1,
#' name = "Poly",
#' D=0.95
#' )
#'
#' final_block=block_rename(2*base_block,c('mu','sigma'))
#'
#' @family {auxiliary functions for structural blocks}
block_rename=function(block,names){
  if(!inherits(block,'dlm_block')){
    stop('Error: The block argument is not a dlm_block object.')
  }
  if(length(names)!=length(block$var_names)){
    stop(paste0('Error: The number of names provided does not match the number of linear predictor in the block. Expected ',length(block$var_names),', got ',length(names),'.'))
  }
  if(length(names)!=length(unique(names))){
    stop(paste0('Error: Repeated names are not allowed.'))
  }

  block$var_names=names
  colnames(block$FF)=names
  return(block)
}
