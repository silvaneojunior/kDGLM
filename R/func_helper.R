#' if.null
#'
#' This function is wrapper for ifelse(is.null(.),.,.)
#'
#' @param vec A vector or matrix.
#' @param val The value to replace NULL with.
#'
#' @export
#' @keywords internal
if.null <- function(vec, val) {
  ifelse(is.null(vec), val, vec)
}

#' if.na
#'
#' This function is wrapper for ifelse(is.na(.),.,.)
#'
#' @param vec A vector or matrix.
#' @param val The value to replace NA with.
#'
#' @export
#' @keywords internal
if.na <- function(vec, val) {
  ifelse(is.na(vec), val, vec)
}


#' if.nan
#'
#' This function is wrapper for ifelse(is.nan(.),.,.)
#'
#' @param vec A vector or matrix.
#' @param val The value to replace NaN with.
#'
#' @export
#' @keywords internal
if.nan <- function(vec, val) {
  ifelse(is.nan(vec), val, vec)
}

#' var_decomp
#'
#' This function receives a covariance matrix S and creates a matrix Q, so that t(Q) %*% Q=S.
#'
#' @param S A covariance matrix
#'
#' @importFrom Rfast cholesky transpose
#' @keywords internal
var_decomp <- function(S) {
  Chol_decomp <- cholesky(S)
  if (prod(if.nan(diag(Chol_decomp), 0)) == 0) {
    svd_decomp <- svd(S)
    d <- sqrt(svd_decomp$d)
    u_t <- transpose(svd_decomp$u)
    return(diag(d) %*% u_t)
  } else {
    return(Chol_decomp)
  }
}

#' ginv
#'
#' This function receives a covariance matrix S and calculates the generalized inverse of S.
#'
#' @param S A covariance matrix
#'
#' @importFrom Rfast cholesky
#' @keywords internal
ginv <- function(S) {
  Chol_decomp <- cholesky(S)
  if (prod(if.nan(diag(Chol_decomp), 0)) < 1e-6) {
    svd_decomp <- svd(S)
    Q_l <- svd_decomp$u
    Q_r <- svd_decomp$v
    D <- svd_decomp$d
    D <- ifelse(D > 1e-6, 1 / D, 0)
    D_mat <- diag(length(D))
    diag(D_mat) <- D
    return(Q_l %*% D_mat %*% transpose(Q_r))
  } else {
    return(chol2inv(Chol_decomp))
  }
}

#' create_G
#'
#' Creates a matrix G such that G %*% S0 %*% t(G)= S1.
#'
#' @param S0 A covariance matrix
#' @param S1 Another covariance matrix
#'
#' @keywords internal
create_G <- function(S0, S1) {
  svd_decomp0 <- svd(S0)
  svd_decomp1 <- svd(S1)
  d0 <- sqrt(svd_decomp0$d)
  d1 <- sqrt(svd_decomp1$d)
  d <- ifelse(d0 > 1e-6, d1 / d0, 0)

  u0 <- transpose(svd_decomp0$u)
  u1 <- svd_decomp1$u
  return(u1 %*% diag(d) %*% u0)
}

#' bdiag
#'
#' Creates a block diagonal matrix with the matrix passed as argument.
#'
#' @param ...  A list of matrices to be used.
#'
#' @return A block diagonal matrix whose diagonal elements are equal to the matrices passed as arguments.
#'
#' @keywords internal
bdiag <- function(...) {
  mats <- list(...)
  ns <- sapply(mats, function(x) {
    if.null(dim(x)[1], 1)
  })
  n <- sum(ns)
  mat_final <- matrix(0, n, n)
  n0 <- 0
  for (mat in mats) {
    n_i <- if.null(dim(mat)[1], 1)
    mat_final[(n0 + 1):(n0 + n_i), (n0 + 1):(n0 + n_i)] <- mat
    n0 <- n0 + n_i
  }
  mat_final
}

#' evaluate_max
#'
#' Auxiliary function to calculate the axis limits and gradation for plots.
#'
#' @param pre_max Numeric: A vector/matrix from which to calculate the axis limits and gradation.
#'
#' @return A list containing the gradation for the axis, the number of ticks in the axis and the maximum value.
#' @keywords internal
evaluate_max <- function(pre_max) {
  if (length(pre_max) == 0 | sum(pre_max**2) < 10**-20) {
    pre_max <- 1
  } else {
    pre_max <- max(pre_max)
  }
  scaled_max <- log10(pre_max)
  category <- scaled_max %% 1
  value <- 10**(floor(log10(max(pre_max))))
  if (category < 0.1) {
    value <- value / 10
  } else {
    if (category < 0.25) {
      value <- value / 5
    } else {
      if (category < 0.5) {
        value <- value / 2
      }
    }
  }
  interval_size <- (pre_max %/% value) + 2
  max_value <- value * interval_size

  return(list(value, interval_size, max_value))
}

#' colQuantile
#'
#' A function that calculates the column-wise quantile of a matrix.
#'
#' @param X Matrix.
#' @param q Numeric: A number between 0 and 1.
#'
#' @importFrom Rfast colnth
#'
#' @export
#' @return Vector: The chosen quantile for each column of X.
#' @keywords internal
colQuantile <- function(X, q) {
  n <- dim(X)[1]
  k <- dim(X)[2]
  min_index <- floor(n * q)
  max_index <- ceiling(n * q)
  (colnth(X, rep(min_index, k)) + colnth(X, rep(max_index, k))) / 2
}

#' rowQuantile
#'
#' A function that calculates the row-wise quantile of a matrix.
#'
#' @param X Matrix.
#' @param q Numeric: A number between 0 and 1.
#'
#' @importFrom Rfast rownth
#'
#' @export
#' @return Vector: The chosen quantile for each row of X.
#' @keywords internal
rowQuantile <- function(X, q) {
  n <- dim(X)[1]
  k <- dim(X)[2]
  min_index <- floor(k * q)
  max_index <- ceiling(k * q)
  (rownth(X, rep(min_index, n)) + rownth(X, rep(max_index, n))) / 2
}

#' f_root
#'
#' Calculates the root of a function given an initial value and a function to calculate it's derivatives.
#'
#' @param f function: A function that receives a vector and return a vector of the same size.
#' @param df function: A function that receives a vector and return the derivatives of f with respect to its arguments (if f returns a vector, it must be a matrix).
#' @param start vector: The initial value to start the algorithm.
#' @param tol numeric: The tolerance for the solution.
#' @param n_max numeric: The maximum number of iterations allowed.
#'
#' @return A list containing:
#' \itemize{
#'    \item root vector: The solution for the system f(x)=0.
#'    \item f.root vector: The function f evaluated at the root.
#'    \item iter numeric: The number of steps taken.
#' }
#'
#' @keywords internal
f_root <- function(f, df, start, tol = 1e-8, n_max = 1000) {
  x_root <- start
  fx <- f(x_root)
  dfx <- df(x_root)
  error <- max(abs(fx))
  count <- 0
  while (error >= tol & count < n_max) {
    count <- count + 1
    change <- solve(dfx, -fx)
    x_root <- x_root + change
    fx <- f(x_root)
    dfx <- df(x_root)
    error <- max(abs(fx))
  }
  if (count >= n_max) {
    warning("Steady state not reached.\n")
  }
  return(list("root" = x_root, "f.root" = fx, "inter." = count))
}

#' check.block.status
#'
#' Checks if a block is defined.
#'
#' @param block A dlm_block object.
#'
#' @return A string ("defined" or "undefined") indicating if all parameters in the block are defined.
#'
#' @import graphics
#'
#' @keywords internal
check.block.status <- function(block) {
  status <- "defined"
  for (param in c("G", "D", "H", "a1", "R1")) {
    if (any(is.character(if.na(block[[param]], 0)))) {
      status <- "undefined"
      break
    }
  }
  if (!all(block$FF_labs %in% c("const", block$pred_names))) {
    status <- "undefined"
  }
  return(status)
}

#' base_ribbon
#'
#' Makes a ribbon plot using R base functions.
#'
#' @param x Vector: A sequence of values for the x-axis.
#' @param ymin Vector: A sequence of values for lower bound of the ribbon.
#' @param ymax Vector: A sequence of values for upper bound of the ribbon.
#' @param ... Extra arguments for the polygon function.
#'
#' @keywords internal
base_ribbon <- function(x, ymin, ymax, ...) {
  l <- length(x)
  polygon(
    x = c(x[1], x, x[l:1]),
    y = c(ymax[1], ymin, ymax[l:1]),
    ...
  )
}
