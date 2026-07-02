################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################
args <- commandArgs(TRUE)
library(R.utils)  
library(stochvol)
library(dplyr)
library(stringr)
library(betategarch)
library(MSGARCH)
library(Rcpp)
library(GAS)
library(rugarch)
source("./DGPs.R")
source("./Utils_GARCH-GAS-SV.R")
sourceCpp("./utils.cpp")
options(scipen = 999)

# Parse arguments if they exist
if (length(args) == 2) {
  dgp <- args[1]
  type <- args[2]
} else {
  dgp <- "FIGARCH-N"
  type <- "US"
}

# Setup
mc <- 1000
n_methods <- 5
n_max <- 10000
n1 <- 3000
n2 <- 5000
n3 <- 10000
garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'norm')
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'std')
figarch_spec_n <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'norm')
figarch_spec_t <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'std')
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
gas_spec_t <- UniGASSpec(Dist = "std",  ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))

if (type == "BR") {
  # data <- logret(read.csv("./Data/precos_diarios_ibrx.csv")[, "PETR4"]) * 100
  garch_params <- c(0.18, 0.09, 0.89)
  figarch_params <- c(0.08, 0.2, 0.5, 0.4)
  dcs_params <- c(0.04, 0.95, 0.05)
  sv_params <- c(1, 0.95, 0.15) # c(0.97, 2.4, 0.16)
  ms_params <- c(0.005, 0.025, 0.95, 0.1, 0.25, 0.70)
  P <-  matrix(c(0.75, 0.30, 0.25, 0.70), 2, 2, byrow = TRUE) 
}
if (type == "US") {
  # data <- readxl::read_xlsx("Data/economatica_nyse_diario.xlsx", skip = 3, col_types = c("text", rep("numeric", 1420)), na = c("-", " ", "NA"))
  # colnames(data) <- stringr::str_replace(colnames(data), "Fechamento\najust p/ prov\nEm moeda orig\n", "")
  # data <- logret(na.omit(data$MMM[3000:10000])) * 100
  garch_params <- c(0.37, 0.14, 0.77)
  figarch_params <- c(0.21, 0.15, 0.3, 0.3)
  dcs_params <- c(0.01, 0.8, 0.1)
  sv_params <- c(1.5, 0.90, 0.20) #c(0.75, 1.5, 0.35) c(1.16, 0.98, 0.12) # c(0.8, 1.7, 0.11) 
  ms_params <- c(0.02, 0.03, 0.93, 0.2, 0.1, 0.80)
  P <- matrix(c(0.95, 0.10, 0.05, 0.90), 2, 2, byrow = TRUE)
}

# Monte Carlo
col_names <- c("returns", "sigma", 
              "garch-n-500", "garch-n-1000", "garch-n-2500",
              "figarch-n-500", "figarch-n-1000", "figarch-n-2500",
              "svb-n-500", "svb-n-1000", "svb-n-2500",
              "ms-n-500", "ms-n-1000", "ms-n-2500",
              "gas-n-500", "gas-n-1000", "gas-n-2500",
              "garch-t-500", "garch-t-1000", "garch-t-2500",
              "figarch-t-500", "figarch-t-1000", "figarch-t-2500",
              "svb-t-500", "svb-t-1000", "svb-t-2500",
              "ms-t-500", "ms-t-1000", "ms-t-2500",
              "gas-t-500", "gas-t-1000", "gas-t-2500")
fore_vols     <- matrix(NA, ncol = 2 + n_methods * 2 * 3, nrow = mc, dimnames = list(NULL, col_names))
fore_vols_o_m <- matrix(NA, ncol = 2 + n_methods * 2 * 3, nrow = mc, dimnames = list(NULL, col_names))
fore_vols_o_w <- matrix(NA, ncol = 2 + n_methods * 2 * 3, nrow = mc, dimnames = list(NULL, col_names))

k <- 0
for (i in 1:mc) {
  print(i)
  is_error <- TRUE
  expr <- NULL
  while (is_error == TRUE) {
    expr <- tryCatch({
        set.seed(i + k)
        data_sim <- switch(dgp,
          "GARCH-N" = garch_sim(n_max + 1, garch_params, "norm"),
          "GARCH-T" = garch_sim(n_max + 1, c(garch_params, 7), "std"),
          "FIGARCH-N" = figarch_sim(n_max + 1, figarch_params, "norm"),
          "FIGARCH-T" = figarch_sim(n_max + 1, c(figarch_params, 7), "std"),
          "DCS-N" = dcs_sim(n_max + 1, dcs_params, "norm"),
          "DCS-T" = dcs_sim(n_max + 1, c(dcs_params, 7), "std"),
          "SV-N" = sv_sim(n_max + 1, sv_params, "norm"),
          "SV-T" = sv_sim(n_max + 1, c(sv_params, 7), "std"),
          "MS-N" = msgarch_sim(n_max + 1, ms_params, "norm", P),
          "MS-T" = msgarch_sim(n_max + 1, c(ms_params, 7), "std", P),
          stop("DGP type not recognized."))
        r_sim <- data_sim$returns
        r_sim_aux <- r_sim
      
        if(any(abs(tail(r_sim, 4)) > 5* sd(r_sim))) stop("Serie mal comportada")
        r_sim_n1 <- tail(r_sim[-c(n_max + 1)], n1)
        r_sim_n2 <- tail(r_sim[-c(n_max + 1)], n2)
        r_sim_n3 <- tail(r_sim[-c(n_max + 1)], n3)
        
        ### No outliers
        garch_n   <- c(as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        figarch_n <- c(as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        gas_n     <- c(sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n1), H = 1)@Forecast$PointForecast[, 2]),
                       sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n2), H = 1)@Forecast$PointForecast[, 2]),
                       sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n3), H = 1)@Forecast$PointForecast[, 2]))
        svb_n     <- c(mean(predict(svsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                       mean(predict(svsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                       mean(predict(svsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
        ms_n      <- c(as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n1), nahead = 1)$vol),
                       as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n2), nahead = 1)$vol),
                       as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n3), nahead = 1)$vol))
        
        garch_t   <- c(as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        figarch_t <- c(as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                       as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
            aux_gas_t_n1 <- gas_fit(gas_spec_t, r_sim_n1)
            nu_n1 <- aux_gas_t_n1@GASDyn$mTheta[3, 1]
            aux_gas_t_n2 <- gas_fit(gas_spec_t, r_sim_n2)
            nu_n2 <- aux_gas_t_n2@GASDyn$mTheta[3, 1]
            aux_gas_t_n3 <- gas_fit(gas_spec_t, r_sim_n3)
            nu_n3 <- aux_gas_t_n3@GASDyn$mTheta[3, 1]
        gas_t     <- c(sqrt(UniGASFor(aux_gas_t_n1, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n1 / (nu_n1 - 2)),
                       sqrt(UniGASFor(aux_gas_t_n2, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n2 / (nu_n2 - 2)),
                       sqrt(UniGASFor(aux_gas_t_n3, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n3 / (nu_n3 - 2)))
        svb_t     <- c(mean(predict(svtsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                       mean(predict(svtsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                       mean(predict(svtsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
        ms_t      <- c(as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n1), nahead = 1)$vol),
                       as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n2), nahead = 1)$vol),
                       as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n3), nahead = 1)$vol))
        
        ### One additive outliers (21 days before the end)
        outlier_position <- n_max + 1 - 22
        r_sim[outlier_position] <- r_sim_aux[outlier_position] + sign(r_sim_aux[outlier_position]) * 5 * sd(r_sim_aux)

        r_sim_n1 <- tail(r_sim[-c(n_max + 1)], n1)
        r_sim_n2 <- tail(r_sim[-c(n_max + 1)], n2)
        r_sim_n3 <- tail(r_sim[-c(n_max + 1)], n3)
        
        garch_n_o_m   <- c(as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        figarch_n_o_m <- c(as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        gas_n_o_m     <- c(sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n1), H = 1)@Forecast$PointForecast[, 2]),
                           sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n2), H = 1)@Forecast$PointForecast[, 2]),
                           sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n3), H = 1)@Forecast$PointForecast[, 2]))
        svb_n_o_m     <- c(mean(predict(svsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
        ms_n_o_m      <- c(as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n1), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n2), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n3), nahead = 1)$vol))
        
        garch_t_o_m   <- c(as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        figarch_t_o_m <- c(as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
            aux_gas_t_n1_o_m <- gas_fit(gas_spec_t, r_sim_n1)
            nu_n1_o_m <- aux_gas_t_n1_o_m@GASDyn$mTheta[3, 1]
            aux_gas_t_n2_o_m <- gas_fit(gas_spec_t, r_sim_n2)
            nu_n2_o_m <- aux_gas_t_n2_o_m@GASDyn$mTheta[3, 1]
            aux_gas_t_n3_o_m <- gas_fit(gas_spec_t, r_sim_n3)
            nu_n3_o_m <- aux_gas_t_n3_o_m@GASDyn$mTheta[3, 1]
        gas_t_o_m     <- c(sqrt(UniGASFor(aux_gas_t_n1_o_m, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n1_o_m / (nu_n1_o_m - 2)),
                           sqrt(UniGASFor(aux_gas_t_n2_o_m, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n2_o_m / (nu_n2_o_m - 2)),
                           sqrt(UniGASFor(aux_gas_t_n3_o_m, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n3_o_m / (nu_n3_o_m - 2)))
        svb_t_o_m     <- c(mean(predict(svtsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svtsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svtsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
        ms_t_o_m      <- c(as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n1), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n2), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n3), nahead = 1)$vol))
        
        ### One additive outliers (5 days before the end)
        outlier_position <- n_max + 1 - 6
        r_sim[outlier_position] <- r_sim_aux[outlier_position] + sign(r_sim_aux[outlier_position]) * 5 * sd(r_sim_aux)
        
        r_sim_n1 <- tail(r_sim[-c(n_max + 1)], n1)
        r_sim_n2 <- tail(r_sim[-c(n_max + 1)], n2)
        r_sim_n3 <- tail(r_sim[-c(n_max + 1)], n3)
        
        garch_n_o_w   <- c(as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_n, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        figarch_n_o_w <- c(as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_n, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        gas_n_o_w     <- c(sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n1), H = 1)@Forecast$PointForecast[, 2]),
                           sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n2), H = 1)@Forecast$PointForecast[, 2]),
                           sqrt(UniGASFor(gas_fit(gas_spec_n, r_sim_n3), H = 1)@Forecast$PointForecast[, 2]))
        svb_n_o_w     <- c(mean(predict(svsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
        ms_n_o_w      <- c(as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n1), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n2), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_n, r_sim_n3), nahead = 1)$vol))
        
        garch_t_o_w   <- c(as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(garch_fit(garch_spec_t, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
        figarch_t_o_w <- c(as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n1), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n2), n.ahead = 1)@forecast$sigmaFor),
                           as.numeric(ugarchforecast(figarch_larger(figarch_spec_t, r_sim_n3), n.ahead = 1)@forecast$sigmaFor))
          aux_gas_t_n1_o_w <- gas_fit(gas_spec_t, r_sim_n1)
          nu_n1_o_w <- aux_gas_t_n1_o_w@GASDyn$mTheta[3, 1]
          aux_gas_t_n2_o_w <- gas_fit(gas_spec_t, r_sim_n2)
          nu_n2_o_w <- aux_gas_t_n2_o_w@GASDyn$mTheta[3, 1]
          aux_gas_t_n3_o_w <- gas_fit(gas_spec_t, r_sim_n3)
          nu_n3_o_w <- aux_gas_t_n3_o_w@GASDyn$mTheta[3, 1]
        gas_t_o_w     <- c(sqrt(UniGASFor(aux_gas_t_n1_o_w, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n1_o_w / (nu_n1_o_w - 2)),
                           sqrt(UniGASFor(aux_gas_t_n2_o_w, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n2_o_w / (nu_n2_o_w - 2)),
                           sqrt(UniGASFor(aux_gas_t_n3_o_w, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu_n3_o_w / (nu_n3_o_w - 2)))
        svb_t_o_w     <- c(mean(predict(svtsample(r_sim_n1, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svtsample(r_sim_n2, quiet = TRUE), steps = 1)$vol[[1]]),
                           mean(predict(svtsample(r_sim_n3, quiet = TRUE), steps = 1)$vol[[1]]))
        ms_t_o_w      <- c(as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n1), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n2), nahead = 1)$vol),
                           as.numeric(predict(msgarch_fit(ms_spec_t, r_sim_n3), nahead = 1)$vol))
        
        fore_vols[i, ] <- c(tail(data_sim$returns, 1), tail(data_sim$volatility, 1), garch_n, figarch_n, svb_n, ms_n, gas_n, garch_t, figarch_t, svb_t, ms_t, gas_t)
        fore_vols_o_m[i, ] <- c(tail(data_sim$returns, 1), tail(data_sim$volatility, 1), garch_n_o_m, figarch_n_o_m, svb_n_o_m, ms_n_o_m, gas_n_o_m, garch_t_o_m, figarch_t_o_m, svb_t_o_m, ms_t_o_m, gas_t_o_m)
        fore_vols_o_w[i, ] <- c(tail(data_sim$returns, 1), tail(data_sim$volatility, 1), garch_n_o_w, figarch_n_o_w, svb_n_o_w, ms_n_o_w, gas_n_o_w, garch_t_o_w, figarch_t_o_w, svb_t_o_w, ms_t_o_w, gas_t_o_w)
        TRUE
      },
      error = function(e) {
        print("Entrou while")
        k <<- k + 10
        return(FALSE)
      })
    if (isTRUE(expr)) {
      is_error <- FALSE
    }
  }
}

write.csv(fore_vols, paste0("./MonteCarlo/volatilities_larger3_", dgp, "_", type, ".csv"))
write.csv(fore_vols_o_m, paste0("./MonteCarlo/volatilities_larger3_outliers_m_", dgp, "_", type, ".csv"))
write.csv(fore_vols_o_w, paste0("./MonteCarlo/volatilities_larger3_outliers_w_", dgp, "_", type, ".csv"))