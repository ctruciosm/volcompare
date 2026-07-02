################################################################################
#####                                DGPs                                 #####
################################################################################


figarch_sim <- function(n, params, distri) {
  # params = (omega, phi, beta, d)
  m <- 1000
  n_burnin <- 5000 + m
  n_tot <- n_burnin + n
  sigma2 <- rep(params[1] / (1 - params[2] - params[3]), n_tot)
  ret <- rep(0, n_tot)
  
  if (distri == "std") {
    epsilon <- rt(n_tot, df = params[5]) * sqrt((params[5] - 2) / params[5])
  } else {
    epsilon <-  rnorm(n_tot)
  }
  
  for (k in 2:n_tot) {
    frac_part <- 0
    for (j in 1:min(k - 1, m)) {
      if (k - 1 - j> 0) {
        frac_part <- frac_part + params[4] * exp(lgamma(j - params[4]) - lgamma(1 - params[4]) - lgamma(j + 1)) * (ret[k - j]^2 - params[2] * ret[k - 1 - j]^2)
      } else {
        frac_part <- frac_part + params[4] * exp(lgamma(j - params[4]) - lgamma(1 - params[4]) - lgamma(j + 1)) * ret[k - j]^2
      }
    }
    sigma2[k] <- params[1] + (params[2] - params[3]) * ret[k - 1]^2 + params[3] * sigma2[k - 1] + frac_part
    ret[k] <- sqrt(sigma2[k]) * epsilon[k]
  }
  return(list(returns = tail(ret, n), volatility = sqrt(tail(sigma2, n)), e = tail(epsilon, n)))
}
#dados <- figarch_sim(10000, c(0.08, 0.2, 0.5, 0.4, 7), "norm")

garch_sim <- function(n, params, distri) { 
  n_burnin <- 5000
  n_tot <- n_burnin + n
  sigma2 <- rep(NA, n_tot)
  ret <- rep(NA, n_tot)
  
  if (distri == "std") {
    epsilon <- rt(n_tot, df = params[4]) * sqrt((params[4] - 2) / params[4])
  } else {
    epsilon <-  rnorm(n_tot)
  }

  sigma2[1] <- params[1] / (1 - params[2] - params[3])
  ret[1] <- sqrt(sigma2[1])*epsilon[1]
  for (i in 2:n_tot) {
    sigma2[i] <- params[1] + params[2] * ret[i - 1]^2 + params[3] *  sigma2[i - 1]
    ret[i] <- sqrt(sigma2[i])*epsilon[i]
  }
  return(list(returns = ret[-c(1:n_burnin)], volatility = sqrt(sigma2[-c(1:n_burnin)]), e = epsilon[-c(1:n_burnin)]))
}
#data <- garch_sim(10000, c(0.01, 0.9, 0.2), "norm")$returns

sv_sim <- function(n, params, distri) {
  # params = [mu, phi, sigma, nu]
  n_burnin <- 5000
  n_tot <- n_burnin + n
  h <- rep(NA, n_tot)
  ret <- rep(NA, n_tot)
  
  eta <- rnorm(n_tot)
  if (distri == "std") {
    epsilon <- rt(n_tot, df = params[4]) * sqrt((params[4] - 2) / params[4])
  } else {
    epsilon <-  rnorm(n_tot)
  }
  
  h[1] <- rnorm(1, params[1], params[3] / sqrt(1 - params[2]^2))
  ret[1] <- exp(h[1]/2)*epsilon[1]
  for (i in 2 : n_tot) {
    h[i] <- params[1] + params[2] * (h[i - 1] - params[1]) + params[3] * eta[i - 1]
    ret[i] <- exp(h[i]/2)*epsilon[i]
  }
  
  vol_one_step_ahead <- exp((params[1] + params[2] * (h[n_tot - 1] - params[1])) / 2) * exp(params[3] ^2 / 8)
  
  return(list(returns = ret[-c(1:n_burnin)], volatility = vol_one_step_ahead, e = epsilon[-c(1:n_burnin)]))
}
  
msgarch_sim <- function(n, params, distri, P) {
  k <- 2
  n_burnin <- 5000
  n_tot <- n_burnin + n 
  h <- matrix(NA, n_tot, k + 1)
  ret <- numeric(n_tot)
  Pt <- numeric(n_tot)
  M <- matrix(NA, 4, 4)
  I4 <- diag(1, 4)
  
  omega <- params[c(1, 4)]
  alpha <- params[c(2, 5)]
  beta  <- params[c(3, 6)]
  if (distri == "std") {
    nu    <- params[7]
    epsilon <- rt(n_tot, df = nu) * sqrt((nu - 2) / nu)
  } else {
    epsilon <-  rnorm(n_tot)
  }
  
  p <- P[1, 1]
  q <- P[2, 2]
  
  M[1, 1] <- P[1, 1] * (alpha[1] + beta[1])
  M[1, 2] <- 0
  M[1, 3] <- P[1, 2] * (alpha[1] + beta[1])
  M[1, 4] <- 0
  M[2, 1] <- P[1, 1] * alpha[2]
  M[2, 2] <- P[1, 1] * beta[2]
  M[2, 3] <- P[1, 2] * alpha[2]
  M[2, 4] <- P[1, 2] * beta[2]
  M[3, 1] <- P[2, 1] * beta[1]
  M[3, 2] <- P[2, 1] * alpha[1]
  M[3, 3] <- P[2, 2] * beta[1]
  M[3, 4] <- P[2, 2] * alpha[1]
  M[4, 1] <- 0;
  M[4, 2] <- P[2, 1] * (alpha[2] + beta[2])
  M[4, 3] <- 0
  M[4, 4] <- P[2, 2] * (alpha[2] + beta[2])
  
  Pt[1] <- (1 - q) / (2 - p - q)       
  pi_inf <- c(Pt[1], 1 - Pt[1])     
  
  h[1, 1:k] <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 2, 4, byrow = TRUE) %*% solve(I4 - M) %*% kronecker(pi_inf, omega)
  h[1, k + 1] <- Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2]
  
  s <- numeric(n_tot)
  s[1] <- 1
  ret[1] <- epsilon[1] * sqrt(h[1, s[1]])
  
  if (distri == "norm") {
    for (i in 2:n_tot) {
      h[i, 1:k] <- omega + alpha * ret[i - 1]^2 + beta * h[i - 1, 1:k]
      Pt[i] <- probability_regime_given_time_n(p, q, sqrt(h[i - 1, 1:k]), ret[i - 1], Pt[i - 1])
      h[i, k + 1] <- Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2]
      s[i] <- sample(1:2, 1, prob = P[, s[i - 1]])
      ret[i] <- epsilon[i] * sqrt(h[i, s[i]])
    }
  } else {
    for (i in 2:n_tot) {
      h[i, 1:k] <- omega + alpha * ret[i - 1]^2 + beta * h[i - 1, 1:k]
      Pt[i] <- probability_regime_given_time_t(p, q, sqrt(h[i - 1, 1:k]), ret[i - 1], Pt[i - 1], 7)
      h[i, k + 1] <- Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2]
      s[i] <- sample(1:2, 1, prob = P[, s[i - 1]])
      ret[i] <- epsilon[i] * sqrt(h[i, s[i]])
    }
  }
  return(list(returns = ret[-c(1:n_burnin)], volatility = sqrt(h[-c(1:n_burnin), k + 1]), e = epsilon[-c(1:n_burnin)], s = s[-c(1:n_burnin)]))
}

dcs_sim <- function(n, params, distri) {
  n_burnin <- 5000
  n_tot <- n_burnin + n
  if (distri == "norm") {
    #Sim <- tegarchSim(n_tot, omega = params[1], phi1 = params[2], phi2 = 0, kappa1 = params[3], kappa2 = 0, kappastar = 0, df = 1200, skew = 1, verbose = TRUE)
    #out <- list(returns = Sim$y[-c(1:n_burnin)], volatility = Sim$stdev[-c(1:n_burnin)], e = Sim$epsilon[-c(1:n_burnin)])
    lambda <- rep(NA, n_tot)
    u <- rep(NA, n_tot)
    ret <- rep(NA, n_tot)
    lambda[1] <- params[1]
    epsilon <-  rnorm(n_tot)
    ret[1] <- exp(lambda[1]) * epsilon[1]
    u[1] <- ret[1]^2 / exp(2 * lambda[1]) - 1
    for (i in 2:n_tot) {
      lambda[i] <- params[1] * (1 - params[2])  + params[2] * lambda[i - 1] + params[3] * u[i - 1]
      ret[i] <- exp(lambda[i]) * epsilon[i]
      u[i] <- ret[i]^2 / exp(2 * lambda[i]) - 1
    }
    out <- list(returns = ret[-c(1:n_burnin)], volatility = exp(lambda[-c(1:n_burnin)]), e = epsilon[-c(1:n_burnin)])
  } else {
    Sim <- tegarchSim(n_tot, omega = params[1], phi1 = params[2], phi2 = 0, kappa1 = params[3], kappa2 = 0, kappastar = 0, df = params[4], skew = 1, verbose = TRUE)
    out <- list(returns = as.numeric(Sim$y[-c(1:n_burnin)]), volatility = as.numeric(Sim$stdev[-c(1:n_burnin)]), e = as.numeric(Sim$epsilon[-c(1:n_burnin)]))
  }
  return(out) 
}




  