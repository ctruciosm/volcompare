################################################################################
#####                   Empirical Application                              #####
################################################################################
library(dplyr)
library(readxl)
library(stringr)
library(rugarch)
library(GAS)
library(stochvol)
library(MSGARCH)
library(future.apply)
source("Utils_GARCH-GAS-SV.R")

# Data
data <- read_excel("./Data/capire_daily_returns.xlsx", skip = 3, col_types = c("date", rep("numeric", 30)), na = c("", "-", NA)) |> 
  filter(Data > "2010-01-01" & Data < "2025-01-01") |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select(where(~ !any(is.na(.x)))) |> 
  rename_with(~ str_remove(.x, "^(?s).*prov\n"), -Data)

ins <- 2500
oos <- nrow(data) - 2500

# Specs
garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
gas_spec_t <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))
figarch_spec_n <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'norm')
figarch_spec_t <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'std')

garch_n_fore_s <- garch_t_fore_s <- figarch_n_fore_s <- figarch_t_fore_s <- gas_n_fore_s <- gas_t_fore_s <- ms_n_fore_s <- ms_t_fore_s <- sv_n_fore_s <- sv_t_fore_s <- svb_n_fore_s <- svb_t_fore_s <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_var1 <- garch_t_fore_var1 <- figarch_n_fore_var1 <- figarch_t_fore_var1 <- gas_n_fore_var1 <- gas_t_fore_var1 <- ms_n_fore_var1 <- ms_t_fore_var1 <- sv_n_fore_var1 <- sv_t_fore_var1 <- svb_n_fore_var1 <- svb_t_fore_var1 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_var2 <- garch_t_fore_var2 <- figarch_n_fore_var2 <- figarch_t_fore_var2 <- gas_n_fore_var2 <- gas_t_fore_var2 <- ms_n_fore_var2 <- ms_t_fore_var2 <- sv_n_fore_var2 <- sv_t_fore_var2 <- svb_n_fore_var2 <- svb_t_fore_var2 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_es1 <- garch_t_fore_es1 <- figarch_n_fore_es1 <- figarch_t_fore_es1 <- gas_n_fore_es1 <- gas_t_fore_es1 <- ms_n_fore_es1 <- ms_t_fore_es1 <- sv_n_fore_es1 <- sv_t_fore_es1 <- svb_n_fore_es1 <- svb_t_fore_es1 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_es2 <- garch_t_fore_es2 <- figarch_n_fore_es2 <- figarch_t_fore_es2 <- gas_n_fore_es2 <- gas_t_fore_es2 <- ms_n_fore_es2 <- ms_t_fore_es2 <- sv_n_fore_es2 <- sv_t_fore_es2 <- svb_n_fore_es2 <- svb_t_fore_es2 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))

#plan(multicore, workers = parallel::detectCores() - 4)
plan(sequential)
for (i in 1:oos) {
  print(i)
  returns <- tail(data[i:(i + 2500 - 1), -1], ins) 
  mu <- apply(returns, 2, mean)
  returns_c <- scale(returns, scale = FALSE)

  for (j in 1:ncol(returns_c)) {
    if(abs(acf(returns_c[, j], plot = FALSE)$acf[2]) > 2/sqrt(nrow(returns_c))) {
      ar_fit <- ar.yw(returns_c[, j], 3, aic = TRUE, se.fit = FALSE)
      mu[j] <- mu[j] + as.numeric(predict(ar_fit, newdata = returns_c[, j], n.ahead = 1, se.fit = FALSE))
      returns_c[, j] <- ar_fit$resid
    }
  }
  
  # Forecasting S2, VaR and ES
  garch_n_fore <- future_apply(returns_c, 2, function(x) {
      fit <- garch_fit(garch_spec_n, na.omit(x)) #ugarchfit(garch_spec_n, na.omit(x), solver = "hybrid")
      alpha <- c(0.01, 0.025)
      e <- fit@fit$residuals / fit@fit$sigma
      sigma_fore <- ugarchforecast(fit, n.ahead = 1)@forecast$sigmaFor[1]
      var_par <- qdist(distribution = "norm", alpha, mu = 0, sigma = sigma_fore)
      es_par <- -sigma_fore * dnorm(qnorm(alpha)) / alpha
      out <- c(sigma_fore, var_par, es_par)  
      return(out)
    }, future.seed = TRUE)
  garch_t_fore <- future_apply(returns_c, 2, function(x) {
    fit <- garch_fit(garch_spec_t, na.omit(x)) #ugarchfit(garch_spec_t, na.omit(x), solver = "hybrid")
    alpha <- c(0.01, 0.025)
    e <- fit@fit$residuals / fit@fit$sigma
    sigma_fore <- ugarchforecast(fit, n.ahead = 1)@forecast$sigmaFor[1]
    nu <- coef(fit)["shape"]
    k <- sqrt(nu / (nu - 2))
    var_par <- sigma_fore * qt(alpha, df = nu) / k
    es_par <- -sigma_fore / k * dt(qt(alpha, df = nu), df = nu) / alpha * (nu + qt(alpha, df = nu)^2) / (nu - 1)
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  figarch_n_fore <- future_apply(returns_c, 2, function(x) {
    fit <- figarch_fit(figarch_spec_n, na.omit(x))
    alpha <- c(0.01, 0.025)
    e <- fit@fit$residuals / fit@fit$sigma
    sigma_fore <- ugarchforecast(fit, n.ahead = 1)@forecast$sigmaFor[1]
    var_par <- qdist(distribution = "norm", alpha, mu = 0, sigma = sigma_fore)
    es_par <- -sigma_fore * dnorm(qnorm(alpha)) / alpha
    out <- c(sigma_fore, var_par, es_par) 
    return(out)
  }, future.seed = TRUE)
  figarch_t_fore <- future_apply(returns_c, 2, function(x) {
    fit <- figarch_fit(figarch_spec_t, na.omit(x))
    alpha <- c(0.01, 0.025)
    e <- fit@fit$residuals / fit@fit$sigma
    sigma_fore <- ugarchforecast(fit, n.ahead = 1)@forecast$sigmaFor[1]
    nu <- coef(fit)["shape"]
    k <- sqrt(nu / (nu - 2))
    var_par <- sigma_fore * qt(alpha, df = nu) / k
    es_par <- -sigma_fore / k * dt(qt(alpha, df = nu), df = nu) / alpha * (nu + qt(alpha, df = nu)^2) / (nu - 1)
    out <- c(sigma_fore, var_par, es_par)
    return(out)
  }, future.seed = TRUE)
  gas_n_fore <- future_apply(returns_c, 2, function(x) {
    fit <- gas_fit(gas_spec_n, na.omit(x))
    alpha <- c(0.01, 0.025)
    sigma <- sqrt(fit@GASDyn$mTheta[2, ])
    e <- fit@Data$vY / sigma[1 : length(fit@Data$vY)]
    sigma_fore <- tail(sigma, 1)
    var_par <- as.numeric(tail(quantile(fit, prob =  alpha), 1))
    es_par <- as.numeric(tail(ES(fit, prob = alpha), 1))
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  gas_t_fore <- future_apply(returns_c, 2, function(x) {
    fit <- gas_fit(gas_spec_t, na.omit(x))
    alpha <- c(0.01, 0.025)
    sigma <- sqrt(fit@GASDyn$mTheta[2, ] * fit@GASDyn$mTheta[3, 1] /(fit@GASDyn$mTheta[3, 1] - 2))
    e <- fit@Data$vY / sigma[1 : length(fit@Data$vY)]
    sigma_fore <- tail(sigma, 1)
    var_par <- as.numeric(tail(quantile(fit, prob =  alpha), 1))
    es_par <- as.numeric(tail(ES(fit, prob = alpha), 1))
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  ms_n_fore <- future_apply(returns_c, 2, function(x) {
    fit <- msgarch_fit(ms_spec_n, na.omit(x))
    alpha <- c(0.01, 0.025)
    e <- fit$data / Volatility(fit)
    sigma_fore <- predict(fit, nahead = 1)$vol
    var_es <- Risk(fit, nahead = 1, alpha = alpha)
    var_par <- var_es$VaR
    es_par <- var_es$ES
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  ms_t_fore <- future_apply(returns_c, 2, function(x) {
    fit <- msgarch_fit(ms_spec_t, na.omit(x))
    alpha <- c(0.01, 0.025)
    e <- fit$data / Volatility(fit)
    sigma_fore <- predict(fit, nahead = 1)$vol
    var_es <- Risk(fit, nahead = 1, alpha = alpha)
    var_par <- var_es$VaR
    es_par <- var_es$ES
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  svb_n_fore <- future_apply(returns_c, 2, function(x) {
    fit <- svsample(as.numeric(na.omit(x)), quiet = TRUE)
    alpha <- c(0.01, 0.025)
    e <- fit$y / as.numeric(colMeans(exp(fit$latent[[1]]/ 2)))
    aux <- predict(fit, steps = 1)
    sigma_fore <- mean(aux$vol[[1]])
    var_par <- quantile(aux$y[[1]], alpha)
    es_par <- sapply(var_par, function(v) mean(aux$y[[1]][aux$y[[1]] < v]))
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  svb_t_fore <- future_apply(returns_c, 2, function(x) {
    fit <- svtsample(as.numeric(na.omit(x)), quiet = TRUE)
    alpha <- c(0.01, 0.025)
    e <- fit$y / as.numeric(colMeans(exp(fit$latent[[1]]/ 2)))
    aux <- predict(fit, steps = 1)
    sigma_fore <- mean(aux$vol[[1]])
    var_par <- quantile(aux$y[[1]], alpha)
    es_par <- sapply(var_par, function(v) mean(aux$y[[1]][aux$y[[1]] < v]))
    out <- c(sigma_fore, var_par, es_par)  
    return(out)
  }, future.seed = TRUE)
  
  # Saving results
  garch_n_fore_s[i, ]    <-  garch_n_fore[1, ]
  garch_t_fore_s[i, ]    <-  garch_t_fore[1, ]
  figarch_n_fore_s[i, ]  <-  figarch_n_fore[1, ]
  figarch_t_fore_s[i, ]  <-  figarch_t_fore[1, ]
  gas_n_fore_s[i, ]      <-  gas_n_fore[1, ]
  gas_t_fore_s[i, ]      <-  gas_t_fore[1, ]
  sv_n_fore_s[i, ]       <-  svb_n_fore[1, ]
  sv_t_fore_s[i, ]       <-  svb_t_fore[1, ]
  ms_n_fore_s[i, ]       <-  ms_n_fore[1, ]
  ms_t_fore_s[i, ]       <-  ms_t_fore[1, ]
  
  garch_n_fore_var1[i, ]    <-  garch_n_fore[2, ] + mu
  garch_t_fore_var1[i, ]    <-  garch_t_fore[2, ] + mu
  figarch_n_fore_var1[i, ]  <-  figarch_n_fore[2, ] + mu
  figarch_t_fore_var1[i, ]  <-  figarch_t_fore[2, ] + mu
  gas_n_fore_var1[i, ]      <-  gas_n_fore[2, ] + mu
  gas_t_fore_var1[i, ]      <-  gas_t_fore[2, ] + mu
  sv_n_fore_var1[i, ]       <-  svb_n_fore[2, ] + mu
  sv_t_fore_var1[i, ]       <-  svb_t_fore[2, ] + mu
  ms_n_fore_var1[i, ]       <-  ms_n_fore[2, ] + mu
  ms_t_fore_var1[i, ]       <-  ms_t_fore[2, ] + mu
  
  garch_n_fore_var2[i, ]    <-  garch_n_fore[3, ] + mu
  garch_t_fore_var2[i, ]    <-  garch_t_fore[3, ] + mu
  figarch_n_fore_var2[i, ]  <-  figarch_n_fore[3, ] + mu
  figarch_t_fore_var2[i, ]  <-  figarch_t_fore[3, ] + mu
  gas_n_fore_var2[i, ]      <-  gas_n_fore[3, ] + mu
  gas_t_fore_var2[i, ]      <-  gas_t_fore[3, ] + mu
  sv_n_fore_var2[i, ]       <-  svb_n_fore[3, ] + mu
  sv_t_fore_var2[i, ]       <-  svb_t_fore[3, ] + mu
  ms_n_fore_var2[i, ]       <-  ms_n_fore[3, ] + mu
  ms_t_fore_var2[i, ]       <-  ms_t_fore[3, ] + mu
  
  garch_n_fore_es1[i, ]    <-  garch_n_fore[4, ] + mu
  garch_t_fore_es1[i, ]    <-  garch_t_fore[4, ] + mu
  figarch_n_fore_es1[i, ]  <-  figarch_n_fore[4, ] + mu
  figarch_t_fore_es1[i, ]  <-  figarch_t_fore[4, ] + mu
  gas_n_fore_es1[i, ]      <-  gas_n_fore[4, ] + mu
  gas_t_fore_es1[i, ]      <-  gas_t_fore[4, ] + mu
  sv_n_fore_es1[i, ]       <-  svb_n_fore[4, ] + mu
  sv_t_fore_es1[i, ]       <-  svb_t_fore[4, ] + mu
  ms_n_fore_es1[i, ]       <-  ms_n_fore[4, ] + mu
  ms_t_fore_es1[i, ]       <-  ms_t_fore[4, ] + mu
  
  garch_n_fore_es2[i, ]    <-  garch_n_fore[5, ] + mu
  garch_t_fore_es2[i, ]    <-  garch_t_fore[5, ] + mu
  figarch_n_fore_es2[i, ]  <-  figarch_n_fore[5, ] + mu
  figarch_t_fore_es2[i, ]  <-  figarch_t_fore[5, ] + mu
  gas_n_fore_es2[i, ]      <-  gas_n_fore[5, ] + mu
  gas_t_fore_es2[i, ]      <-  gas_t_fore[5, ] + mu
  sv_n_fore_es2[i, ]       <-  svb_n_fore[5, ] + mu
  sv_t_fore_es2[i, ]       <-  svb_t_fore[5, ] + mu
  ms_n_fore_es2[i, ]       <-  ms_n_fore[5, ] + mu
  ms_t_fore_es2[i, ]       <-  ms_t_fore[5, ] + mu
}


write.csv(garch_n_fore_s, "Empirical_Application/garch_n_fore_s_2500.csv")
write.csv(garch_t_fore_s, "Empirical_Application/garch_t_fore_s_2500.csv")
write.csv(figarch_n_fore_s, "Empirical_Application/figarch_n_fore_s_2500.csv")
write.csv(figarch_t_fore_s, "Empirical_Application/figarch_t_fore_s_2500.csv")
write.csv(gas_n_fore_s, "Empirical_Application/gas_n_fore_s_2500.csv")
write.csv(gas_t_fore_s, "Empirical_Application/gas_t_fore_s_2500.csv")
write.csv(ms_n_fore_s, "Empirical_Application/ms_n_fore_s_2500.csv")
write.csv(ms_t_fore_s, "Empirical_Application/ms_t_fore_s_2500.csv")
write.csv(sv_n_fore_s, "Empirical_Application/sv_n_fore_s_2500.csv")
write.csv(sv_t_fore_s, "Empirical_Application/sv_t_fore_s_2500.csv")

write.csv(garch_n_fore_var1, "Empirical_Application/garch_n_fore_var1_2500.csv")
write.csv(garch_t_fore_var1, "Empirical_Application/garch_t_fore_var1_2500.csv")
write.csv(figarch_n_fore_var1, "Empirical_Application/figarch_n_fore_var1_2500.csv")
write.csv(figarch_t_fore_var1, "Empirical_Application/figarch_t_fore_var1_2500.csv")
write.csv(gas_n_fore_var1, "Empirical_Application/gas_n_fore_var1_2500.csv")
write.csv(gas_t_fore_var1, "Empirical_Application/gas_t_fore_var1_2500.csv")
write.csv(ms_n_fore_var1, "Empirical_Application/ms_n_fore_var1_2500.csv")
write.csv(ms_t_fore_var1, "Empirical_Application/ms_t_fore_var1_2500.csv")
write.csv(sv_n_fore_var1, "Empirical_Application/sv_n_fore_var1_2500.csv")
write.csv(sv_t_fore_var1, "Empirical_Application/sv_t_fore_var1_2500.csv")
  
write.csv(garch_n_fore_var2, "Empirical_Application/garch_n_fore_var2_2500.csv")
write.csv(garch_t_fore_var2, "Empirical_Application/garch_t_fore_var2_2500.csv")
write.csv(figarch_n_fore_var2, "Empirical_Application/figarch_n_fore_var2_2500.csv")
write.csv(figarch_t_fore_var2, "Empirical_Application/figarch_t_fore_var2_2500.csv")
write.csv(gas_n_fore_var2, "Empirical_Application/gas_n_fore_var2_2500.csv")
write.csv(gas_t_fore_var2, "Empirical_Application/gas_t_fore_var2_2500.csv")
write.csv(ms_n_fore_var2, "Empirical_Application/ms_n_fore_var2_2500.csv")
write.csv(ms_t_fore_var2, "Empirical_Application/ms_t_fore_var2_2500.csv")
write.csv(sv_n_fore_var2, "Empirical_Application/sv_n_fore_var2_2500.csv")
write.csv(sv_t_fore_var2, "Empirical_Application/sv_t_fore_var2_2500.csv")

write.csv(garch_n_fore_es1, "Empirical_Application/garch_n_fore_es1_2500.csv")
write.csv(garch_t_fore_es1, "Empirical_Application/garch_t_fore_es1_2500.csv")
write.csv(figarch_n_fore_es1, "Empirical_Application/figarch_n_fore_es1_2500.csv")
write.csv(figarch_t_fore_es1, "Empirical_Application/figarch_t_fore_es1_2500.csv")
write.csv(gas_n_fore_es1, "Empirical_Application/gas_n_fore_es1_2500.csv")
write.csv(gas_t_fore_es1, "Empirical_Application/gas_t_fore_es1_2500.csv")
write.csv(ms_n_fore_es1, "Empirical_Application/ms_n_fore_es1_2500.csv")
write.csv(ms_t_fore_es1, "Empirical_Application/ms_t_fore_es1_2500.csv")
write.csv(sv_n_fore_es1, "Empirical_Application/sv_n_fore_es1_2500.csv")
write.csv(sv_t_fore_es1, "Empirical_Application/sv_t_fore_es1_2500.csv")

write.csv(garch_n_fore_es2, "Empirical_Application/garch_n_fore_es2_2500.csv")
write.csv(garch_t_fore_es2, "Empirical_Application/garch_t_fore_es2_2500.csv")
write.csv(figarch_n_fore_es2, "Empirical_Application/figarch_n_fore_es2_2500.csv")
write.csv(figarch_t_fore_es2, "Empirical_Application/figarch_t_fore_es2_2500.csv")
write.csv(gas_n_fore_es2, "Empirical_Application/gas_n_fore_es2_2500.csv")
write.csv(gas_t_fore_es2, "Empirical_Application/gas_t_fore_es2_2500.csv")
write.csv(ms_n_fore_es2, "Empirical_Application/ms_n_fore_es2_2500.csv")
write.csv(ms_t_fore_es2, "Empirical_Application/ms_t_fore_es2_2500.csv")
write.csv(sv_n_fore_es2, "Empirical_Application/sv_n_fore_es2_2500.csv")
write.csv(sv_t_fore_es2, "Empirical_Application/sv_t_fore_es2_2500.csv")

write.csv(data[(1 + 2500):nrow(data), ], "Empirical_Application/data_oos_2500.csv")