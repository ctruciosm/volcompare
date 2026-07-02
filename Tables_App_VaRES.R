################################################################################
#####                  Tables Empirical Application                        #####
################################################################################
library(dplyr)
library(lubridate)
library(modelconf)
library(wrMisc)
library(readxl)
library(kableExtra)
library(GAS)
library(esback)
library(rugarch)
library(esreg)
library(Rcpp)
source("Utils_GARCH-GAS-SV.R")
sourceCpp("utils.cpp")


# Import Data
datas_oos <- read.csv("./Empirical_Application/data_oos_2500.csv")[, -1]

s_garch_n <- read.csv("./Empirical_Application/garch_n_fore_s_2500.csv")[, -1]
s_garch_t <- read.csv("./Empirical_Application/garch_t_fore_s_2500.csv")[, -1]
s_figarch_n <- read.csv("./Empirical_Application/figarch_n_fore_s_2500.csv")[, -1]
s_figarch_t <- read.csv("./Empirical_Application/figarch_t_fore_s_2500.csv")[, -1]
s_gas_n <- read.csv("./Empirical_Application/gas_n_fore_s_2500.csv")[, -1]
s_gas_t <- read.csv("./Empirical_Application/gas_t_fore_s_2500.csv")[, -1]
s_ms_n <- read.csv("./Empirical_Application/ms_n_fore_s_2500.csv")[, -1]
s_ms_t <- read.csv("./Empirical_Application/ms_t_fore_s_2500.csv")[, -1]
s_sv_n <- read.csv("./Empirical_Application/sv_n_fore_s_2500.csv")[, -1]
s_sv_t <- read.csv("./Empirical_Application/sv_t_fore_s_2500.csv")[, -1]


var1_garch_n <- read.csv("./Empirical_Application/garch_n_fore_var1_2500.csv")[, -1]
var1_garch_t <- read.csv("./Empirical_Application/garch_t_fore_var1_2500.csv")[, -1]
var1_figarch_n <- read.csv("./Empirical_Application/figarch_n_fore_var1_2500.csv")[, -1]
var1_figarch_t <- read.csv("./Empirical_Application/figarch_t_fore_var1_2500.csv")[, -1]
var1_gas_n <- read.csv("./Empirical_Application/gas_n_fore_var1_2500.csv")[, -1]
var1_gas_t <- read.csv("./Empirical_Application/gas_t_fore_var1_2500.csv")[, -1]
var1_ms_n <- read.csv("./Empirical_Application/ms_n_fore_var1_2500.csv")[, -1]
var1_ms_t <- read.csv("./Empirical_Application/ms_t_fore_var1_2500.csv")[, -1]
var1_sv_n <- read.csv("./Empirical_Application/sv_n_fore_var1_2500.csv")[, -1]
var1_sv_t <- read.csv("./Empirical_Application/sv_t_fore_var1_2500.csv")[, -1]

es1_garch_n <- read.csv("./Empirical_Application/garch_n_fore_es1_2500.csv")[, -1]
es1_garch_t <- read.csv("./Empirical_Application/garch_t_fore_es1_2500.csv")[, -1]
es1_figarch_n <- read.csv("./Empirical_Application/figarch_n_fore_es1_2500.csv")[, -1]
es1_figarch_t <- read.csv("./Empirical_Application/figarch_t_fore_es1_2500.csv")[, -1]
es1_gas_n <- read.csv("./Empirical_Application/gas_n_fore_es1_2500.csv")[, -1]
es1_gas_t <- read.csv("./Empirical_Application/gas_t_fore_es1_2500.csv")[, -1]
es1_ms_n <- read.csv("./Empirical_Application/ms_n_fore_es1_2500.csv")[, -1]
es1_ms_t <- read.csv("./Empirical_Application/ms_t_fore_es1_2500.csv")[, -1]
es1_sv_n <- read.csv("./Empirical_Application/sv_n_fore_es1_2500.csv")[, -1]
es1_sv_t <- read.csv("./Empirical_Application/sv_t_fore_es1_2500.csv")[, -1]

var2_garch_n <- read.csv("./Empirical_Application/garch_n_fore_var2_2500.csv")[, -1]
var2_garch_t <- read.csv("./Empirical_Application/garch_t_fore_var2_2500.csv")[, -1]
var2_figarch_n <- read.csv("./Empirical_Application/figarch_n_fore_var2_2500.csv")[, -1]
var2_figarch_t <- read.csv("./Empirical_Application/figarch_t_fore_var2_2500.csv")[, -1]
var2_gas_n <- read.csv("./Empirical_Application/gas_n_fore_var2_2500.csv")[, -1]
var2_gas_t <- read.csv("./Empirical_Application/gas_t_fore_var2_2500.csv")[, -1]
var2_ms_n <- read.csv("./Empirical_Application/ms_n_fore_var2_2500.csv")[, -1]
var2_ms_t <- read.csv("./Empirical_Application/ms_t_fore_var2_2500.csv")[, -1]
var2_sv_n <- read.csv("./Empirical_Application/sv_n_fore_var2_2500.csv")[, -1]
var2_sv_t <- read.csv("./Empirical_Application/sv_t_fore_var2_2500.csv")[, -1]

es2_garch_n <- read.csv("./Empirical_Application/garch_n_fore_es2_2500.csv")[, -1]
es2_garch_t <- read.csv("./Empirical_Application/garch_t_fore_es2_2500.csv")[, -1]
es2_figarch_n <- read.csv("./Empirical_Application/figarch_n_fore_es2_2500.csv")[, -1]
es2_figarch_t <- read.csv("./Empirical_Application/figarch_t_fore_es2_2500.csv")[, -1]
es2_gas_n <- read.csv("./Empirical_Application/gas_n_fore_es2_2500.csv")[, -1]
es2_gas_t <- read.csv("./Empirical_Application/gas_t_fore_es2_2500.csv")[, -1]
es2_ms_n <- read.csv("./Empirical_Application/ms_n_fore_es2_2500.csv")[, -1]
es2_ms_t <- read.csv("./Empirical_Application/ms_t_fore_es2_2500.csv")[, -1]
es2_sv_n <- read.csv("./Empirical_Application/sv_n_fore_es2_2500.csv")[, -1]
es2_sv_t <- read.csv("./Empirical_Application/sv_t_fore_es2_2500.csv")[, -1]



# Back Testing
K <- ncol(es2_garch_n)
a1 <- 0.010
a2 <- 0.025
out <- NULL
out_1 <- NULL
out_2 <- NULL
for (i in 1:K) { 
  print(i)
  set.seed(12345)
  r_oos <- datas_oos[,1 + i]
  
  hits1 <- round(c(mean(r_oos < var1_garch_n[, i]), mean(r_oos < var1_garch_t[, i]), 
            mean(r_oos < var1_figarch_n[, i]), mean(r_oos < var1_figarch_t[, i]), 
            mean(r_oos < var1_gas_n[, i]), mean(r_oos < var1_gas_t[, i]), 
            mean(r_oos < var1_ms_n[, i]), mean(r_oos < var1_ms_t[, i]), 
            mean(r_oos < var1_sv_n[, i]), mean(r_oos < var1_sv_t[, i]))*100, 2)
  tests1 <- c(calibration_tests(r_oos, var1_garch_n[, i], es1_garch_n[, i], s_garch_n[, i], alpha = a1),
              calibration_tests(r_oos, var1_garch_t[, i], es1_garch_t[, i], s_garch_t[, i], alpha = a1),
              calibration_tests(r_oos, var1_figarch_n[, i], es1_figarch_n[, i], s_figarch_n[, i], alpha = a1),
              calibration_tests(r_oos, var1_figarch_t[, i], es1_figarch_t[, i], s_figarch_t[, i], alpha = a1),
              calibration_tests(r_oos, var1_gas_n[, i], es1_gas_n[, i], s_gas_n[, i], alpha = a1),
              calibration_tests(r_oos, var1_gas_t[, i], es1_gas_t[, i], s_gas_t[, i], alpha = a1),
              calibration_tests(r_oos, var1_ms_n[, i], es1_ms_n[, i], s_ms_n[, i], alpha = a1),
              calibration_tests(r_oos, var1_ms_t[, i], es1_ms_t[, i], s_ms_t[, i], alpha = a1),
              calibration_tests(r_oos, var1_sv_n[, i], es1_sv_n[, i], s_sv_n[, i], alpha = a1),
              calibration_tests(r_oos, var1_sv_t[, i], es1_sv_t[, i], s_sv_t[, i], alpha = a1))
  
  scores1 <- mcs_scoring_functions(r_oos, 
                                  var = cbind(var1_garch_n[, i], var1_garch_t[, i], var1_figarch_n[, i], var1_figarch_t[, i], var1_gas_n[, i], var1_gas_t[, i], var1_ms_n[, i], var1_ms_t[, i], var1_sv_n[, i], var1_sv_t[, i]),
                                  es = cbind(es1_garch_n[, i], es1_garch_t[, i], es1_figarch_n[, i], es1_figarch_t[, i], es1_gas_n[, i], es1_gas_t[, i], es1_ms_n[, i], es1_ms_t[, i], es1_sv_n[, i], es1_sv_t[, i]),
                                  alpha = a1)
  
  hits2 <- round(c(mean(r_oos < var2_garch_n[, i]), mean(r_oos < var2_garch_t[, i]), 
                   mean(r_oos < var2_figarch_n[, i]), mean(r_oos < var2_figarch_t[, i]), 
                   mean(r_oos < var2_gas_n[, i]), mean(r_oos < var2_gas_t[, i]), 
                   mean(r_oos < var2_ms_n[, i]), mean(r_oos < var2_ms_t[, i]), 
                   mean(r_oos < var2_sv_n[, i]), mean(r_oos < var2_sv_t[, i]))*100, 2)
  tests2 <- round(c(calibration_tests(r_oos, var2_garch_n[, i], es2_garch_n[, i], s_garch_n[, i], alpha = a2),
              calibration_tests(r_oos, var2_garch_t[, i], es2_garch_t[, i], s_garch_t[, i], alpha = a2),
              calibration_tests(r_oos, var2_figarch_n[, i], es2_figarch_n[, i], s_figarch_n[, i], alpha = a2),
              calibration_tests(r_oos, var2_figarch_t[, i], es2_figarch_t[, i], s_figarch_t[, i], alpha = a2),
              calibration_tests(r_oos, var2_gas_n[, i], es2_gas_n[, i], s_gas_n[, i], alpha = a2),
              calibration_tests(r_oos, var2_gas_t[, i], es2_gas_t[, i], s_gas_t[, i], alpha = a2),
              calibration_tests(r_oos, var2_ms_n[, i], es2_ms_n[, i], s_ms_n[, i], alpha = a2),
              calibration_tests(r_oos, var2_ms_t[, i], es2_ms_t[, i], s_ms_t[, i], alpha = a2),
              calibration_tests(r_oos, var2_sv_n[, i], es2_sv_n[, i], s_sv_n[, i], alpha = a2),
              calibration_tests(r_oos, var2_sv_t[, i], es2_sv_t[, i], s_sv_t[, i], alpha = a2)), 1)

  
  scores2 <- round(mcs_scoring_functions(r_oos, 
                                   var = cbind(var2_garch_n[, i], var2_garch_t[, i], var2_figarch_n[, i], var2_figarch_t[, i], var2_gas_n[, i], var2_gas_t[, i], var2_ms_n[, i], var2_ms_t[, i], var2_sv_n[, i], var2_sv_t[, i]),
                                   es = cbind(es2_garch_n[, i], es2_garch_t[, i], es2_figarch_n[, i], es2_figarch_t[, i], es2_gas_n[, i], es2_gas_t[, i], es2_ms_n[, i], es2_ms_t[, i], es2_sv_n[, i], es2_sv_t[, i]),
                                   alpha = a2), 1)
  
  aux <- rbind(c(hits1, hits2), c(tests1, tests2), c(scores1, scores2))
  out <- rbind(out, aux)
  out_1 <- rbind(out_1, if_else(tests1 * scores1 == 16, 1, 0))
  out_2 <- rbind(out_2, if_else(tests2 * scores2 == 16, 1, 0))
}

group_names <- colnames(es1_garch_n)
n_per_group <- 3
models_names <- c("GARCH-N", "GARCH-T", "FIGARCH-N", "FIGARCH-T", "GAS-N", "GAS-T", "MS-N", "MS-T", "SV-N", "SV-T")

table_latex <- out |> 
  kbl(format = "latex", booktabs = TRUE, escape = FALSE, col.names = rep(models_names, 2)) |> 
  add_header_above(c(" " = 1, "1\\%" = 10, "2.5\\%" = 10)) 

start <- seq(1, nrow(out), by = n_per_group)
end   <- pmin(start + n_per_group - 1, nrow(out))

for (i in seq_along(start)) {
  table_latex <- table_latex |>
    pack_rows(group_names[i], start[i], end[i])
}

table_latex |> save_kable(file = "App_VaRES_2500.tex", keep_tex = TRUE)

