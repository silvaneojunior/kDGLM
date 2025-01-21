#' if.null
#'
#' This function is wrapper for ifelse(is.null(.),.,.)
#'
#' @param vec A vector or matrix.
#' @param val The value to replace NULL with.
#'
#' @return A vector or matrix with the same dimensions as the input, where any NULL values have been replaced by the specified val argument.
#'
#' @export
#' @keywords internal
if.null <- function(vec, val) {
  ifelse(is.null(vec), val, vec)
}

#' if.na
#'
#' This function is wrapper for ifelse(is.na(vec),vec,val)
#'
#' @param vec A vector or matrix.
#' @param val The value to replace NA with.
#'
#' @return A vector or matrix with the same dimensions as the input, where any NA values have been replaced by the specified val argument.
#'
#' @export
#' @keywords internal
if.na <- function(vec, val) {
  ifelse(is.na(vec), val, vec)
}


#' if.nan
#'
#' This function is wrapper for ifelse(is.nan(vec),vec,val)
#'
#' @param vec A vector or matrix.
#' @param val The value to replace NaN with.
#'
#' @return A vector or matrix with the same dimensions as the input, where any NaN values have been replaced by the specified val argument.
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
#' @importFrom Rfast cholesky transpose rowAny
#' @keywords internal
var_decomp <- function(S) {
  # S=as.matrix(S)
  # diag(S)=diag(S)+1e-8
  Chol.decomp <- cholesky(S)
  flags.index <- is.nan(Chol.decomp) | is.infinite(Chol.decomp)
  flags <- rowAny(flags.index)
  while (any(flags)) {
    sub.chol <- Chol.decomp[!flags, flags, drop = FALSE]
    placeholder <- S[flags, flags, drop = FALSE] - t(sub.chol) %*% sub.chol
    # diag(placeholder)=abs(diag(placeholder))
    # flags.neg=diag(placeholder)<1e-6
    # placeholder[flags.neg,]=0
    # order=order(diag(placeholder))
    # flags2 <- diag(placeholder) > 1e-6
    flags2 <- diag(placeholder) > 0
    sub.flags <- flags
    sub.flags[flags] <- flags2
    if (any(sub.flags)) {
      Chol.decomp[sub.flags, sub.flags] <- cholesky(placeholder[flags2, flags2])
    }
    flags.index <- is.nan(Chol.decomp) | is.infinite(Chol.decomp)
    Chol.decomp[flags.index] <- 0
    flags <- rowAny(flags.index)
  }
  return(Chol.decomp)
}

# var_decomp <- function(S) {
#   Chol.decomp <- cholesky(S)
#   if (prod(diag(Chol.decomp), na.rm = TRUE) == 0 || any(is.na(Chol.decomp))) {
#     svd.decomp <- svd(S)
#     d <- sqrt(svd.decomp$d)
#     u.t <- transpose(svd.decomp$u)
#     return(diag(d, nrow = length(d)) %*% u.t)
#   } else {
#     return(Chol.decomp)
#   }
# }


# var_decomp3 <- function(S) {
#   Chol.decomp <- cholesky(S)
#   flags=is.nan(diag(Chol.decomp))
#   while(any(flags)){
#     Chol.decomp[flags,flags]=cholesky(S[flags,flags])
#     Chol.decomp[which(flags)[1]-1,]=0
#     flags=is.nan(diag(Chol.decomp))
#   }
#   return(Chol.decomp)
# }

# var_decomp2 <- function(S) {
#   Chol.decomp <- cholesky(S)
#   if (prod(diag(Chol.decomp), na.rm = TRUE) == 0 || any(is.na(Chol.decomp))) {
#     svd.decomp <- svd(S)
#     d <- sqrt(svd.decomp$d)
#     u.t <- transpose(svd.decomp$u)
#     return(diag(d, nrow = length(d)) %*% u.t)
#   } else {
#     return(Chol.decomp)
#   }
# }

#' ginv
#'
#' This function receives a covariance matrix S and calculates the generalized inverse of S.
#'
#' @param S A covariance matrix
#'
#' @importFrom Rfast cholesky
#' @keywords internal
ginv <- function(S) {
  S <- as.matrix(S)
  Chol.decomp <- var_decomp(S)
  flags <- diag(Chol.decomp) > 1e-8
  if (!all(flags)) {
    inv <- S * 0
    inv[flags, flags] <- chol2inv(Chol.decomp[flags, flags])
  } else {
    inv <- chol2inv(Chol.decomp)
  }
  return(inv)
}
# ginv <- function(S) {
#   Chol.decomp <- cholesky(S)
#   if (prod(diag(Chol.decomp), na.rm = TRUE) < 1e-12 || any(is.na(Chol.decomp))) {
#     svd.decomp <- svd(S)
#     Q.l <- svd.decomp$u
#     Q.r <- svd.decomp$v
#     D <- svd.decomp$d
#     D <- ifelse(D > 1e-12, 1 / D, 0)
#     D.mat <- diag(length(D))
#     diag(D.mat) <- D
#     return(Q.l %*% D.mat %*% transpose(Q.r))
#   } else {
#     return(chol2inv(Chol.decomp))
#   }
# }



#' dmvnorm
#'
#' Calculates the log density of a multivariate normal distribution with mean mu and covariance matrix Sigma.
#'
#' @param x Vector: The value from to which calculate the density.
#' @param mu Vector: The mean vector
#' @param Sigma Matrix: The Covariance matrix.
#'
#' @importFrom stats dnorm
#'
#' @keywords internal
dmvnorm <- function(x, mu, Sigma) {
  Sigma <- as.matrix(Sigma)
  index <- order(diag(Sigma), decreasing = TRUE)
  x <- x[index]
  mu <- mu[index]
  Sigma <- Sigma[index, index, drop = FALSE]
  diag(Sigma) <- diag(Sigma) + 1e-8


  Chol.decomp <- var_decomp(Sigma)
  flags.valid <- diag(Chol.decomp) > 0
  inv.chol.Sigma <- Chol.decomp * 0

  inv.chol.Sigma[flags.valid, flags.valid] <-
    backsolve(Chol.decomp[flags.valid, flags.valid], diag(sum(flags.valid)))
  diag.chol.inv <- diag(inv.chol.Sigma)

  norm.x <- transpose(inv.chol.Sigma) %*% (x - mu)

  # Chol.decomp <- cholesky(Sigma)
  # if (prod(diag(Chol.decomp), na.rm = TRUE) < 1e-12 || any(is.na(Chol.decomp))) {
  #   svd.decomp <- svd(Sigma)
  #   Q.l <- svd.decomp$u
  #   Q.r <- svd.decomp$v
  #   D <- svd.decomp$d
  #   D <- ifelse(D > 1e-12, 1 / sqrt(D), 0)
  #   D.mat <- diag(length(D))
  #   diag(D.mat) <- D
  #   inv.chol.Sigma <- (Q.l %*% D.mat %*% transpose(Q.r))
  #   diag.chol.inv <- D
  # } else {
  #   inv.chol.Sigma <- (backsolve(Chol.decomp, diag(length(mu))))
  #   diag.chol.inv <- diag(inv.chol.Sigma)
  # }
  #
  # flags.valid <- diag.chol.inv > 1e-12
  # norm.x <- transpose(inv.chol.Sigma) %*% (x - mu)

  sum(dnorm(norm.x[flags.valid], log = TRUE)) +
    sum(log(abs(diag.chol.inv[flags.valid])))
}


#' rmvnorm
#'
#' Obtains a sample from a multivariate normal distribution.
#'
#' @param n integer: The sample size.
#' @param mu numeric: The mean vector
#' @param Sigma matrix: The Covariance matrix.
#'
#' @importFrom Rfast matrnorm
#' @importFrom stats runif
#'
#' @keywords internal
rmvnorm <- function(n, mu, Sigma,
                    norm.x = matrnorm(k, n, seed = round(runif(1) * 1e15))) {
  k <- length(mu)
  chol.Sigma <- var_decomp(Sigma)
  transpose(chol.Sigma) %*% norm.x + c(mu)
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
  svd.decomp0 <- solve(var_decomp(S0))
  svd.decomp1 <- var_decomp(S1)
  return(transpose(svd.decomp0 %*% svd.decomp1))
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
  mat.final <- matrix(0, n, n)
  n.0 <- 0
  for (mat in mats) {
    n.i <- if.null(dim(mat)[1], 1)
    mat.final[(n.0 + 1):(n.0 + n.i), (n.0 + 1):(n.0 + n.i)] <- mat
    n.0 <- n.0 + n.i
  }
  mat.final
}

#' evaluate_max
#'
#' Auxiliary function to calculate the axis limits and gradation for plots.
#'
#' @param pre.max numeric: A vector/matrix from which to calculate the axis limits and gradation.
#'
#' @return A list containing the gradation for the axis, the number of ticks in the axis and the maximum value.
#' @keywords internal
evaluate_max <- function(pre.max) {
  if (length(pre.max) == 0 || sum(pre.max**2) < 10**-20) {
    pre.max <- 1
  } else {
    pre.max <- max(pre.max)
  }
  scaled.max <- log10(pre.max)
  category <- scaled.max %% 1
  value <- 10**(floor(log10(max(pre.max))))
  if (category < 0.1) {
    value <- value / 10
  } else if (category < 0.25) {
    value <- value / 5
  } else if (category < 0.5) {
    value <- value / 2
  }
  interval.size <- (pre.max %/% value) + 2
  max.value <- value * interval.size

  return(list(value, interval.size, max.value))
}

#' colQuantile
#'
#' A function that calculates the column-wise quantile of a matrix.
#'
#' @param X matrix.
#' @param q numeric: A number between 0 and 1.
#'
#' @importFrom Rfast colnth
#'
#' @export
#' @return numeric: The chosen quantile for each column of X.
#' @keywords internal
colQuantile <- function(X, q) {
  n <- dim(X)[1]
  k <- dim(X)[2]
  min.index <- floor(n * q)
  max.index <- ceiling(n * q)
  (colnth(X, rep(min.index, k)) + colnth(X, rep(max.index, k))) / 2
}

#' rowQuantile
#'
#' A function that calculates the row-wise quantile of a matrix.
#'
#' @param X matrix.
#' @param q numeric: A number between 0 and 1.
#'
#' @importFrom Rfast rownth
#'
#' @export
#' @return numeric: The chosen quantile for each row of X.
#' @keywords internal
rowQuantile <- function(X, q) {
  n <- dim(X)[1]
  k <- dim(X)[2]
  min.index <- floor(k * q)
  max.index <- ceiling(k * q)
  (rownth(X, rep(min.index, n)) + rownth(X, rep(max.index, n))) / 2
}

#' f_root
#'
#' Calculates the root of a function given an initial value and a function to calculate its derivatives.
#'
#' @param f function: A function that receives a vector and return a vector of the same size.
#' @param df function: A function that receives a vector and return the derivatives of f with respect to its arguments (if f returns a vector, it must be a matrix).
#' @param start vector: The initial value to start the algorithm.
#' @param tol numeric: The tolerance for the solution.
#' @param n.max numeric: The maximum number of iterations allowed.
#'
#' @return A list containing:
#' \itemize{
#'    \item root vector: The solution for the system f(x)=0.
#'    \item f.root vector: The function f evaluated at the root.
#'    \item iter numeric: The number of steps taken.
#' }
#'
#' @keywords internal
f_root <- function(f, df, start, tol = 1e-8, n.max = 1000) {
  k <- length(start)
  S <- tol * diag(k)
  x.root <- start
  fx <- f(x.root)
  dfx <- df(x.root) + S
  error <- max(abs(fx))
  count <- 0
  while (error >= tol && count < n.max) {
    count <- count + 1
    # print('#####################')
    # print(det(dfx))
    # print(fx)
    change <- solve(dfx, -fx)
    # print(change)
    x.root <- x.root + change
    fx <- f(x.root)
    dfx <- df(x.root) + S
    error <- max(abs(fx))
  }
  if (count >= n.max) {
    warning("Steady state not reached.\n")
  }
  return(list("root" = x.root, "f.root" = fx, "inter." = count))
}

#' f_joint_root
#'
#' Calculates the root of a function given an initial value and a function to calculate its derivatives.
#'
#' @param f function: A function that receives a vector and return a vector of the same size and a matrix representing its derivatives.
#' @param start vector: The initial value to start the algorithm.
#' @param tol numeric: The tolerance for the solution.
#' @param n.max numeric: The maximum number of iterations allowed.
#'
#' @return A list containing:
#' \itemize{
#'    \item root vector: The solution for the system f(x)=0.
#'    \item f.root vector: The function f evaluated at the root.
#'    \item iter numeric: The number of steps taken.
#' }
#'
#' @keywords internal
f_joint_root <- function(f, start, tol = 1e-8, n.max = 1000) {
  k <- length(start)
  S <- 0.001 * diag(k)
  x.root <- start
  fx.joint <- f(x.root)
  fx <- fx.joint[[1]]
  dfx <- fx.joint[[2]] + S
  error <- max(abs(fx))
  count <- 0
  while (error >= tol && count < n.max) {
    count <- count + 1
    change <- solve(dfx, -fx)
    x.root <- x.root + change
    fx.joint <- f(x.root)
    fx <- fx.joint[[1]]
    dfx <- fx.joint[[2]] + S
    error <- max(abs(fx))
  }
  if (count >= n.max) {
    warning("Steady state not reached.\n")
  }
  return(list("root" = x.root, "f.root" = fx, "inter." = count))
}

#' check.block.status
#'
#' Checks if a block is defined.
#'
#' @param block A dlm_block object.
#'
#' @return A character ("defined" or "undefined") indicating if all parameters in the block are defined.
#'
#' @import graphics
#'
#' @keywords internal
check.block.status <- function(block) {
  status <- "defined"
  for (param in c("G", "D", "H", "a1", "R1", "scale")) {
    if (any(is.character(if.na(block[[param]], 0)))) {
      status <- "undefined"
      break
    }
  }
  if (!all(block$FF.labs %in% c("const", "Covariate", block$pred.names))) {
    status <- "undefined"
  }
  return(status)
}

#' base_ribbon
#'
#' Makes a ribbon plot using R base functions.
#'
#' @param x numeric: A sequence of values for the x-axis.
#' @param ymin numeric: A sequence of values for lower bound of the ribbon.
#' @param ymax numeric: A sequence of values for upper bound of the ribbon.
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

#' lcm
#'
#' Calculates the least common multiple of a set of integer. Internal use only.
#'
#' @param x numeric: A sequence of integers.
#'
#' @return The least common multiple.
#'
#' @keywords internal
lcm <- function(x) {
  if (any(x != round(x))) {
    return(1)
  }
  if (length(x) == 2) {
    y <- prod(x)
    vals1 <- seq.int(x[2], y, x[2])
    return(min(vals1[(vals1 %% x[1]) == 0]))
  } else {
    return(lcm(c(x[1], lcm(x[-1]))))
  }
}

# na.conditions

all_na <- function(data, offset) {
  all(is.na(data) | is.na(offset) | offset == 0)
}
any_na <- function(data, offset) {
  any(is.na(data) | is.na(offset) | offset == 0)
}
