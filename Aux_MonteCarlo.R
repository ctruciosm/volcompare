################################################################################
#####              Auxiliary Monte Carlo Experiments                       #####
################################################################################
library(R.utils)
library(xtable)
library(betategarch)
library(stochvol)
library(stochvolTMB)
library(dplyr)
library(stringr)
library(MSGARCH)
library(Rcpp)
library(GAS)
library(forecast)
library(kableExtra)
library(modelconf)
source("./DGPs.R")
source("./Utils_GARCH-GAS-SV.R")
sourceCpp("./utils.cpp")
options(scipen = 999)

# Setup
mc <- 500
n_methods <- 7
n_max <- 2500
n1 <- 500
n2 <- 1000
n3 <- 2500
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))
gas_spec_t <- UniGASSpec(Dist = "std",  ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))

#if (type == "BR") {
#  dcs_params <- c(0.04, 0.95, 0.05)
#  sv_params <- c(1, 0.95, 0.15) # c(0.97, 2.4, 0.16)
#  ms_params <- c(0.005, 0.025, 0.95, 0.1, 0.25, 0.70)
#  P <-  matrix(c(0.75, 0.30, 0.25, 0.70), 2, 2, byrow = TRUE) 
#}
#if (type == "US") {
dcs_params <- c(0.01, 0.8, 0.1)
sv_params <- c(1.5, 0.90, 0.20) #c(0.75, 1.5, 0.35) c(1.16, 0.98, 0.12) # c(0.8, 1.7, 0.11) 
ms_params <- c(0.02, 0.03, 0.93, 0.2, 0.1, 0.80)
P <- matrix(c(0.95, 0.10, 0.05, 0.90), 2, 2, byrow = TRUE)
#}


# GAS
fore_vols_gas  <- matrix(NA, ncol = 7, nrow = mc, dimnames = list(NULL, c("Sigma", "GAS-500", "GAS-1000", "GAS-2500", "DCS-500", "DCS-1000", "DCS-2500")))
for (i in 1:mc) {
  set.seed(i)
  print(i)
  data_sim <- dcs_sim(n_max + 1, c(dcs_params, 7), "std")
  is_error <- TRUE
  expr <- NULL
  while (is_error == TRUE) {
    r_sim <- data_sim$returns
    r_sim_aux <- r_sim
    expr <- tryCatch({
      r_sim_n1 <- tail(r_sim[-c(n_max + 1)], n1)
      r_sim_n2 <- tail(r_sim[-c(n_max + 1)], n2)
      r_sim_n3 <- tail(r_sim[-c(n_max + 1)], n3)
      # Package betategarch
      dcs_t <- c(as.numeric(predict(tegarch(r_sim_n1, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev),
                 as.numeric(predict(tegarch(r_sim_n2, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev),
                 as.numeric(predict(tegarch(r_sim_n3, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev))
      # Package GAS
      aux_gas_t_n1 <- gas_fit(gas_spec_t, r_sim_n1)
        nu_n1 <- aux_gas_t_n1@GASDyn$mTheta[3, 1]
      aux_gas_t_n2 <- gas_fit(gas_spec_t, r_sim_n2)
        nu_n2 <- aux_gas_t_n2@GASDyn$mTheta[3, 1]
      aux_gas_t_n3 <- gas_fit(gas_spec_t, r_sim_n3)
        nu_n3 <- aux_gas_t_n3@GASDyn$mTheta[3, 1]
      gas_t     <- c(sqrt(UniGASFor(aux_gas_t_n1, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n1 / (nu_n1 - 2)),
                     sqrt(UniGASFor(aux_gas_t_n2, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n2 / (nu_n2 - 2)),
                     sqrt(UniGASFor(aux_gas_t_n3, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n3 / (nu_n3 - 2)))
      TRUE
    },
      error = function(e) {
        data_sim <- dcs_sim(n_max + 1, c(dcs_params, 7), "std")
        FALSE
      })
    if (isTRUE(expr)) {
      is_error <- FALSE
    }
  }
  fore_vols_gas[i, ] <- c(tail(data_sim$volatility, 1), gas_t, dcs_t)
}
write.csv(fore_vols_gas, "auxiliary_fore_vols_gas.csv")

# Diebold-Mariano test for GAS vs betategarch
e1 <- fore_vols_gas[, "Sigma"] - fore_vols_gas[, "GAS-500"]
e2 <- fore_vols_gas[, "Sigma"] - fore_vols_gas[, "DCS-500"]
dm_gas_500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_gas[, "Sigma"] - fore_vols_gas[, "GAS-1000"]
e2 <- fore_vols_gas[, "Sigma"] - fore_vols_gas[, "DCS-1000"]
dm_gas_1000 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_gas[, "Sigma"] - fore_vols_gas[, "GAS-2500"]
e2 <- fore_vols_gas[, "Sigma"] - fore_vols_gas[, "DCS-2500"]
dm_gas_2500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)


# MSGARCH
fore_vols_ms  <- matrix(NA, ncol = 7, nrow = mc, dimnames = list(NULL, c("Sigma", "MSB-500", "MSB-1000", "MSB-2500", "MS-500", "MS-1000", "MS-2500")))
for (i in 1:mc) {
  set.seed(i)
  print(i)
  data_sim <- msgarch_sim(n_max + 1, c(ms_params, 7), "std", P)
  is_error <- TRUE
  expr <- NULL
  while (is_error == TRUE) {
    r_sim <- data_sim$returns
    r_sim_aux <- r_sim
    expr <- tryCatch({
      r_sim_n1 <- tail(r_sim[-c(n_max + 1)], n1)
      r_sim_n2 <- tail(r_sim[-c(n_max + 1)], n2)
      r_sim_n3 <- tail(r_sim[-c(n_max + 1)], n3)
      # Package MSGARCH MCMC
      msb <- c(as.numeric(predict(FitMCMC(ms_spec_t, r_sim_n1), nahead = 1)$vol),
        as.numeric(predict(FitMCMC(ms_spec_t, r_sim_n2), nahead = 1)$vol),
        as.numeric(predict(FitMCMC(ms_spec_t, r_sim_n3), nahead = 1)$vol))
      # Package MSGARCH ML
      ms <- c(as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n1), nahead = 1)$vol),
        as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n2), nahead = 1)$vol),
        as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n3), nahead = 1)$vol))
      TRUE
    },
      error = function(e) {
        data_sim <- msgarch_sim(n_max + 1, c(ms_params, 7), "std", P)
        FALSE
      })
    if (isTRUE(expr)) {
      is_error <- FALSE
    }
  }
  fore_vols_ms[i, ] <- c(tail(data_sim$volatility, 1), msb, ms)
}
write.csv(fore_vols_ms, "auxiliary_fore_vols_ms.csv")

# Diebold-Mariano test for GAS vs betategarch
e1 <- fore_vols_ms[, "Sigma"] - fore_vols_ms[, "MSB-500"]
e2 <- fore_vols_ms[, "Sigma"] - fore_vols_ms[, "MS-500"]
dm_ms_500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_ms[, "Sigma"] - fore_vols_ms[, "MSB-1000"]
e2 <- fore_vols_ms[, "Sigma"] - fore_vols_ms[, "MS-1000"]
dm_ms_1000 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_ms[, "Sigma"] - fore_vols_ms[, "MSB-2500"]
e2 <- fore_vols_ms[, "Sigma"] - fore_vols_ms[, "MS-2500"]
dm_ms_2500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)


# SV
fore_vols_sv  <- matrix(NA, ncol = 7, nrow = mc, dimnames = list(NULL, c("Sigma", "SV-500", "SV-1000", "SV-2500", "TMB-500", "TMB-1000", "TMB-2500")))
for (i in 1:mc) {
  set.seed(i)
  print(i)
  data_sim <- sv_sim(n_max + 1, c(sv_params, 7), "std")
  is_error <- TRUE
  expr <- NULL
  while (is_error == TRUE) {
    r_sim <- data_sim$returns
    r_sim_aux <- r_sim
    expr <- tryCatch({
      r_sim_n1 <- tail(r_sim[-c(n_max + 1)], n1)
      r_sim_n2 <- tail(r_sim[-c(n_max + 1)], n2)
      r_sim_n3 <- tail(r_sim[-c(n_max + 1)], n3)
      # Package stochvol
      svb_t     <- c(mean(predict(svtsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                     mean(predict(svtsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                     mean(predict(svtsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
      # Package stochvolTMB
      aux_sv_t_n1 <- predict(estimate_parameters_t(r_sim_n1), steps = 1)$h_exp
      aux_sv_t_n2 <- predict(estimate_parameters_t(r_sim_n2), steps = 1)$h_exp
      aux_sv_t_n3 <- predict(estimate_parameters_t(r_sim_n3), steps = 1)$h_exp
      sv_t      <- c(mean(aux_sv_t_n1), mean(aux_sv_t_n2), mean(aux_sv_t_n3))
      TRUE
    },
      error = function(e) {
        data_sim <- sv_sim(n_max + 1, c(sv_params, 7), "std")
        FALSE
      })
    if (isTRUE(expr)) {
      is_error <- FALSE
    }
  }
  fore_vols_sv[i, ] <- c(tail(data_sim$volatility, 1), svb_t, sv_t)
}
write.csv(fore_vols_sv, "auxiliary_fore_vols_sv.csv")

# Diebold-Mariano test for GAS vs betategarch
e1 <- fore_vols_sv[, "Sigma"] - fore_vols_sv[, "SV-500"]
e2 <- fore_vols_sv[, "Sigma"] - fore_vols_sv[, "TMB-500"]
dm_sv_500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_sv[, "Sigma"] - fore_vols_sv[, "SV-1000"]
e2 <- fore_vols_sv[, "Sigma"] - fore_vols_sv[, "TMB-1000"]
dm_sv_1000 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_sv[, "Sigma"] - fore_vols_sv[, "SV-2500"]
e2 <- fore_vols_sv[, "Sigma"] - fore_vols_sv[, "TMB-2500"]
dm_sv_2500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)


# Tables
table_latex <- rbind(round(c(mean(loss_mse(fore_vols_gas[, "Sigma"], fore_vols_gas[, "GAS-500"])),
                           mean(loss_mse(fore_vols_gas[, "Sigma"], fore_vols_gas[, "DCS-500"])),
                           as.numeric(dm_gas_500$p.value),
                           mean(loss_mse(fore_vols_gas[, "Sigma"], fore_vols_gas[, "GAS-1000"])),
                           mean(loss_mse(fore_vols_gas[, "Sigma"], fore_vols_gas[, "DCS-1000"])),
                           as.numeric(dm_gas_1000$p.value),
                           mean(loss_mse(fore_vols_gas[, "Sigma"], fore_vols_gas[, "GAS-2500"])),
                           mean(loss_mse(fore_vols_gas[, "Sigma"], fore_vols_gas[, "DCS-2500"])),
                           as.numeric(dm_gas_2500$p.value)), 4),
                     round(c(mean(loss_mse(fore_vols_ms[, "Sigma"], fore_vols_ms[, "MSB-500"])),
                           mean(loss_mse(fore_vols_ms[, "Sigma"], fore_vols_ms[, "MS-500"])),
                           as.numeric(dm_ms_500$p.value),
                           mean(loss_mse(fore_vols_ms[, "Sigma"], fore_vols_ms[, "MSB-1000"])),
                           mean(loss_mse(fore_vols_ms[, "Sigma"], fore_vols_ms[, "MS-1000"])),
                           as.numeric(dm_ms_1000$p.value),
                           mean(loss_mse(fore_vols_ms[, "Sigma"], fore_vols_ms[, "MSB-2500"])),
                           mean(loss_mse(fore_vols_ms[, "Sigma"], fore_vols_ms[, "MS-2500"])),
                           as.numeric(dm_ms_2500$p.value)), 4),
                     round(c(mean(loss_mse(fore_vols_sv[, "Sigma"], fore_vols_sv[, "SV-500"])),
                           mean(loss_mse(fore_vols_sv[, "Sigma"], fore_vols_sv[, "TMB-500"])),
                           as.numeric(dm_sv_500$p.value),
                           mean(loss_mse(fore_vols_sv[, "Sigma"], fore_vols_sv[, "SV-1000"])),
                           mean(loss_mse(fore_vols_sv[, "Sigma"], fore_vols_sv[, "TMB-1000"])),
                           as.numeric(dm_sv_1000$p.value),
                           mean(loss_mse(fore_vols_sv[, "Sigma"], fore_vols_sv[, "SV-2500"])),
                           mean(loss_mse(fore_vols_sv[, "Sigma"], fore_vols_sv[, "TMB-2500"])),
                           as.numeric(dm_sv_2500$p.value)), 4))

M <- as.matrix(table_latex)
rownames(M) <- c("GAS", "MSGARCH", "SV")
colnames(M) <- c("Classic-500", "Alternative-500", "p-value-500",
                 "Classic-1000", "Alternative-1000", "p-value-1000",
                 "Classic-2500", "Alternative-2500", "p-value-2500")
print(xtable(M, digits = 4), file = "auxiliary_monte_carlo.tex", include.rownames = TRUE, include.colnames = TRUE)
                   



fore_vols_gas1 <- read.csv("auxiliary_fore_vols_gas.csv")[, -1]
fore_vols_ms1  <- read.csv("auxiliary_fore_vols_ms.csv")[, -1]
fore_vols_sv1  <- read.csv("auxiliary_fore_vols_sv.csv")[, -1]

e1 <- fore_vols_gas1[, "Sigma"] - fore_vols_gas1[, "GAS.500"]
e2 <- fore_vols_gas1[, "Sigma"] - fore_vols_gas1[, "DCS.500"]
dm_gas_500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_gas1[, "Sigma"] - fore_vols_gas1[, "GAS.1000"]
e2 <- fore_vols_gas1[, "Sigma"] - fore_vols_gas1[, "DCS.1000"]
dm_gas_1000 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_gas1[, "Sigma"] - fore_vols_gas1[, "GAS.2500"]
e2 <- fore_vols_gas1[, "Sigma"] - fore_vols_gas1[, "DCS.2500"]
dm_gas_2500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)


e1 <- fore_vols_ms1[, "Sigma"] - fore_vols_ms1[, "MSB.500"]
e2 <- fore_vols_ms1[, "Sigma"] - fore_vols_ms1[, "MS.500"]
dm_ms_500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_ms1[, "Sigma"] - fore_vols_ms1[, "MSB.1000"]
e2 <- fore_vols_ms1[, "Sigma"] - fore_vols_ms1[, "MS.1000"]
dm_ms_1000 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_ms1[, "Sigma"] - fore_vols_ms1[, "MSB.2500"]
e2 <- fore_vols_ms1[, "Sigma"] - fore_vols_ms1[, "MS.2500"]
dm_ms_2500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)

e1 <- fore_vols_sv1[, "Sigma"] - fore_vols_sv1[, "SV.500"]
e2 <- fore_vols_sv1[, "Sigma"] - fore_vols_sv1[, "TMB.500"]
dm_sv_500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_sv1[, "Sigma"] - fore_vols_sv1[, "SV.1000"]
e2 <- fore_vols_sv1[, "Sigma"] - fore_vols_sv1[, "TMB.1000"]
dm_sv_1000 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)
e1 <- fore_vols_sv1[, "Sigma"] - fore_vols_sv1[, "SV.2500"]
e2 <- fore_vols_sv1[, "Sigma"] - fore_vols_sv1[, "TMB.2500"]
dm_sv_2500 <- dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2)

table_latex <- rbind(round(c(mean(loss_mse(fore_vols_gas1[, "Sigma"], fore_vols_gas1[, "GAS.500"])),
  mean(loss_mse(fore_vols_gas1[, "Sigma"], fore_vols_gas1[, "DCS.500"])),
  as.numeric(dm_gas_500$p.value),
  mean(loss_mse(fore_vols_gas1[, "Sigma"], fore_vols_gas1[, "GAS.1000"])),
  mean(loss_mse(fore_vols_gas1[, "Sigma"], fore_vols_gas1[, "DCS.1000"])),
  as.numeric(dm_gas_1000$p.value),
  mean(loss_mse(fore_vols_gas1[, "Sigma"], fore_vols_gas1[, "GAS.2500"])),
  mean(loss_mse(fore_vols_gas1[, "Sigma"], fore_vols_gas1[, "DCS.2500"])),
  as.numeric(dm_gas_2500$p.value)), 4),
  round(c(mean(loss_mse(fore_vols_ms1[, "Sigma"], fore_vols_ms1[, "MSB.500"])),
    mean(loss_mse(fore_vols_ms1[, "Sigma"], fore_vols_ms1[, "MS.500"])),
    as.numeric(dm_ms_500$p.value),
    mean(loss_mse(fore_vols_ms1[, "Sigma"], fore_vols_ms1[, "MSB.1000"])),
    mean(loss_mse(fore_vols_ms1[, "Sigma"], fore_vols_ms1[, "MS.1000"])),
    as.numeric(dm_ms_1000$p.value),
    mean(loss_mse(fore_vols_ms1[, "Sigma"], fore_vols_ms1[, "MSB.2500"])),
    mean(loss_mse(fore_vols_ms1[, "Sigma"], fore_vols_ms1[, "MS.2500"])),
    as.numeric(dm_ms_2500$p.value)), 4),
  round(c(mean(loss_mse(fore_vols_sv1[, "Sigma"], fore_vols_sv1[, "SV.500"])),
    mean(loss_mse(fore_vols_sv1[, "Sigma"], fore_vols_sv1[, "TMB.500"])),
    as.numeric(dm_sv_500$p.value),
    mean(loss_mse(fore_vols_sv1[, "Sigma"], fore_vols_sv1[, "SV.1000"])),
    mean(loss_mse(fore_vols_sv1[, "Sigma"], fore_vols_sv1[, "TMB.1000"])),
    as.numeric(dm_sv_1000$p.value),
    mean(loss_mse(fore_vols_sv1[, "Sigma"], fore_vols_sv1[, "SV.2500"])),
    mean(loss_mse(fore_vols_sv1[, "Sigma"], fore_vols_sv1[, "TMB.2500"])),
    as.numeric(dm_sv_2500$p.value)), 4))

M <- as.matrix(table_latex)
rownames(M) <- c("GAS", "MSGARCH", "SV")
colnames(M) <- c("Classic-500", "Alternative-500", "p-value-500",
  "Classic-1000", "Alternative-1000", "p-value-1000",
  "Classic-2500", "Alternative-2500", "p-value-2500")
print(xtable(M, digits = 4), file = "auxiliary_monte_carlo.tex", include.rownames = TRUE, include.colnames = TRUE)


# Tables for larger sample sizes

garch_n <- read.csv("./MonteCarlo/volatilities_GARCH-N_BR.csv")[, -1]^2
garch_t <- read.csv("./MonteCarlo/volatilities_GARCH-T_BR.csv")[, -1]^2
figarch_n <- read.csv("./Aux_MonteCarlo/volatilities_larger_FIGARCH-N_BR.csv")[, -1]^2
figarch_t <- read.csv("./Aux_MonteCarlo/volatilities_larger_FIGARCH-T_BR.csv")[, -1]^2
dcs_n <- read.csv("./MonteCarlo/volatilities_DCS-N_BR.csv")[, -1]^2
dcs_t <- read.csv("./MonteCarlo/volatilities_DCS-T_BR.csv")[, -1]^2
sv_n <- read.csv("./MonteCarlo/volatilities_SV-N_BR.csv")[, -1]^2
sv_t <- read.csv("./MonteCarlo/volatilities_SV-T_BR.csv")[, -1]^2
ms_n <- read.csv("./Aux_MonteCarlo/volatilities_larger_MS-N_BR.csv")[, -1]^2
ms_t <- read.csv("./Aux_MonteCarlo/volatilities_larger_MS-T_BR.csv")[, -1]^2

garch_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_GARCH-N_BR.csv")[, -1]^2
garch_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_GARCH-T_BR.csv")[, -1]^2
figarch_n_out_m <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_m_FIGARCH-N_BR.csv")[, -1]^2
figarch_t_out_m <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_m_FIGARCH-T_BR.csv")[, -1]^2
dcs_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_DCS-N_BR.csv")[, -1]^2
dcs_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_DCS-T_BR.csv")[, -1]^2
sv_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_SV-N_BR.csv")[, -1]^2
sv_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_SV-T_BR.csv")[, -1]^2
ms_n_out_m <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_m_MS-N_BR.csv")[, -1]^2
ms_t_out_m <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_m_MS-T_BR.csv")[, -1]^2

garch_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_GARCH-N_BR.csv")[, -1]^2
garch_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_GARCH-T_BR.csv")[, -1]^2
figarch_n_out_w <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_w_FIGARCH-N_BR.csv")[, -1]^2
figarch_t_out_w <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_w_FIGARCH-T_BR.csv")[, -1]^2
dcs_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_DCS-N_BR.csv")[, -1]^2
dcs_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_DCS-T_BR.csv")[, -1]^2
sv_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_SV-N_BR.csv")[, -1]^2
sv_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_SV-T_BR.csv")[, -1]^2
ms_n_out_w <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_w_MS-N_BR.csv")[, -1]^2
ms_t_out_w <- read.csv("./Aux_MonteCarlo/volatilities_larger_outliers_w_MS-T_BR.csv")[, -1]^2



# Tables MSE, QLIKE, MAE
set.seed(1234)
for (l in c("MSE", "QLIKE")) {
  for (i in c(0)) {
    table_latex(loss_func = l, w = i) |> save_kable(file = paste0("tex_tables/", l, "_", i, "_BR_75_Large.tex"), keep_tex = TRUE)
  }
}

