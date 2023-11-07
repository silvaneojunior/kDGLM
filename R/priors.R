#' Zero sum prior
#'
#' Set the prior of a structural block to be such that the latent variables sum zero.
#' The covariance matrix of the evolution and the drift parameter are also altered to guarantee that the zero sum condition will always hold.
#'
#' @param block dlm_block: The structural block.
#' @param var.index Integer: The index of the variables from which to set the prior.
#'
#' @return A dlm_block object with the desired prior.
#'
#' @export
#' @examples
#' @details
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' @family {auxiliary functions for structural blocks}
#'
#' @references
#'    \insertAllCited{}
zero_sum_prior <- function(block, var.index = 1:block$n) {
  transf <- matrix(-1 / block$n, block$n, block$n)
  diag(transf) <- 1 + diag(transf)
  block$a1 <- block$a1 - mean(block$a1)
  block$R1 <- transf %*% block$R1 %*% transf
  for (i in 1:block$t) {
    d <- min(block$D[, , i])
    if (min(block$D[, , i]) > 0 & min(block$D[, , i]) != max(block$D[, , i])) {
      d <- block$D[, , i]
      d <- mean(d[d != 0])
      warning("Not all latent states have the same discount factor. All values will be set to the average value.")
    }
    block$D[, , i] <- d
    block$h[, i] <- block$h[, i] - mean(block$h[, i])
    block$H[, , i] <- transf %*% block$H[, , i] %*% transf
  }
  return(block)
}

#' CAR prior
#'
#' Set the prior of a structural block as a Conditional Autoregressive (CAR) prior.
#'
#' @param block dlm_block: The structural block.
#' @param adj.matrix Matrix: The adjacency matrix.
#' @param tau Numeric: The tau parameter for the CAR model (see references).
#' @param rho Numeric: The rho parameter for the CAR model (see references).
#' @return A dlm_block object with the desired prior.
#'
#' @importFrom Rfast is.symmetric
#' @export
#' @examples
#' @details
#'
#' For a revision of the CAR prior, see \insertCite{AlexCar;textual}{kDGLM}.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' @family {auxiliary functions for structural blocks}
#'
#' @references
#'    \insertAllCited{}
CAR_prior <- function(block, adj.matrix, scale, rho, var.index = 1:block$n) {
  if (dim(adj.matrix)[1] != length(var.index)) {
    stop("Error: Number of regions must be equal to the number of variables.")
  }
  k <- dim(adj.matrix)[1]
  adj.matrix <- matrix(adj.matrix, k, k)
  if (!is.symmetric(adj.matrix)) {
    stop("Error: adj.matrix is not symmetric.")
  }
  D.mat <- diag(rowSums(adj.matrix))
  R <- (D.mat - adj.matrix)
  R1 <- ginv(((1 - rho) * diag(k) + rho * R))
  if (is.character(scale)) {
    block$scale <- scale
  } else {
    R1 <- R1 * scale
  }
  block$R1 <- R1

  return(block)
}
