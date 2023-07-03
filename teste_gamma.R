
system_full_gamma <- function(x, parms) {
  n <- exp(x) # exp(x[1])
  k <- n # exp(x[2])
  tau <- (n * parms$Hq1 + 1) / parms$Hq2
  theta <- n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1)

  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # print((parms$Hq3 + parms$Hq4))
  if (a <= 5) {
    # Densidade marginal aproximada de alpha (uso opcional).
    # f_densi=function(x){dgamma(a,b))}
    # c_val=1
    # Densidade marginal exata de phi.
    f_densi_raw <- function(x) {
      exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
    }
    lim_sup <- Inf
    c_val <- cubintegrate(f_densi_raw, 0, lim_sup, nVec = 200)$integral
    f_densi <- function(x) {
      f_densi_raw(x) / c_val
    }
    # print('a')
    f <- function(x) {
      (x * digamma(x * n + 1) - x * log(x) - x * log(tau)) * f_densi(x)
    }
    Hp3 <- cubintegrate(f, 0, lim_sup, nVec = 200)$integral

    # print('b')
    f <- function(x) {
      (x * log(x) - lgamma(x)) * f_densi(x)
    }
    Hp4 <- cubintegrate(f, 0, lim_sup, nVec = 200)$integral

    # print('c')
    # f <- function(x) {
    #   (x * digamma(x * n + 1)  - lgamma(x)- x*log(tau)) * f_densi(x)
    # }
    # Hp5 <- cubintegrate(f, 0, Inf, nVec = 200)$integral
    # print('sd')
    Hp5 <- Hp3 + Hp4


    # f <- function(x) {
    #   (lgamma(x)-x*digamma(n*x+1)) * f_densi(x)
    # }
    # Hp5 <- integrate(f, 0, Inf)$value+parms$Hq1*log(tau)
  } else {
    c_val <- 1
    # Hp3 <- log(tau / n) * a / b - 1 / n + b / (12 * (n**2) * (a - 1))
    # Hp4 <- a / b + 0.5 * (digamma(a) - log(b)) - b / (12 * (a - 1)) - 11 / 12
    # Hp5=Hp3+Hp4
    Hp5 <- parms$Hq1 * (log(tau / n) - 1) - 0.5 * digamma((n + 1) / 2) + 0.5 * log(n * log(tau / n) - theta) + (log(tau / n) - theta / n) / 6 + 11 / 12 - 1 / (2 * n)
  }

  f_all <- c(
    (Hp5 - (parms$Hq3 + parms$Hq4))
  )
  # print(f_all)
  return(f_all)
}


convert_FGamma_Normal <- function(ft, Qt, parms) {
  s <- exp(ft[2, ] - 1)
  f1 <- ft[1, ]
  f2 <- ft[2, ] - log(s)
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  Hq1 <- exp(f1 + q1 / 2)
  Hq2 <- exp(f1 - f2 + (q1 + q2 - 2 * q12) / 2)
  Hq3 <- -(f2 - f1*q12/q1) * Hq1+(q12/q1)*(f1+q1)*exp(((f1+q1)**2)/(2*q1)-((1**2)/(2*q1)))

  Hq4 <- cubintegrate(function(x) {
    (x * log(x) - lgamma(x)) * dlnorm(x, f1, sqrt(q1))
  }, 0, Inf, nVec = 200)$integral

  parms <- list(
    "Hq1" = Hq1,
    "Hq2" = Hq2,
    "Hq3" = Hq3,
    "Hq4" = Hq4
  )

  ss1 <- multiroot(f = system_full_gamma, start = c(0), parms = parms, maxiter = 2000, atol = 10**-20)

  x <- as.numeric(ss1$root)
  n <- exp(x) # exp(x[1])

  k <- n # exp(x[2])

  # Calculando tau e theta dado n e k
  tau <- ((n * parms$Hq1 + 1) / parms$Hq2)
  # tau=exp(x[3])
  theta <- (n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1))
  tau <- tau * s
  theta <- theta + n * log(s)
  return(list("n" = n, "k" = k, "tau" = tau, "theta" = theta))


}
