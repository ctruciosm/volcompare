################################################################################
#####                       Tables Monte Carlo                             #####
################################################################################
library(dplyr)
library(ggplot2)
library(tidyr)
library(modelconf)
library(wrMisc)
library(xtable)
library(kableExtra)
source("Utils_GARCH-GAS-SV.R")
options(scipen = 999)
# Importing

garch_n <- read.csv("./MonteCarlo/volatilities_GARCH-N_BR.csv")[, -1]^2
garch_t <- read.csv("./MonteCarlo/volatilities_GARCH-T_BR.csv")[, -1]^2
figarch_n <- read.csv("./Aux_MonteCarlo/volatilities_larger2_FIGARCH-N_BR.csv")[, -1]^2
figarch_t <- read.csv("./Aux_MonteCarlo/volatilities_larger2_FIGARCH-T_BR.csv")[, -1]^2
dcs_n <- read.csv("./MonteCarlo/volatilities_DCS-N_BR.csv")[, -1]^2
dcs_t <- read.csv("./MonteCarlo/volatilities_DCS-T_BR.csv")[, -1]^2
sv_n <- read.csv("./MonteCarlo/volatilities_SV-N_BR.csv")[, -1]^2
sv_t <- read.csv("./MonteCarlo/volatilities_SV-T_BR.csv")[, -1]^2
ms_n <- read.csv("./MonteCarlo/volatilities_MS-N_BR.csv")[, -1]^2
ms_t <- read.csv("./MonteCarlo/volatilities_MS-T_BR.csv")[, -1]^2

garch_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_GARCH-N_BR.csv")[, -1]^2
garch_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_GARCH-T_BR.csv")[, -1]^2
figarch_n_out_m <- read.csv("./Aux_MonteCarlo/volatilities_larger2_outliers_m_FIGARCH-N_BR.csv")[, -1]^2
figarch_t_out_m <- read.csv("./Aux_MonteCarlo/volatilities_larger2_outliers_m_FIGARCH-T_BR.csv")[, -1]^2
dcs_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_DCS-N_BR.csv")[, -1]^2
dcs_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_DCS-T_BR.csv")[, -1]^2
sv_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_SV-N_BR.csv")[, -1]^2
sv_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_SV-T_BR.csv")[, -1]^2
ms_n_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_MS-N_BR.csv")[, -1]^2
ms_t_out_m <- read.csv("./MonteCarlo/volatilities_outliers_m_MS-T_BR.csv")[, -1]^2

garch_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_GARCH-N_BR.csv")[, -1]^2
garch_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_GARCH-T_BR.csv")[, -1]^2
figarch_n_out_w <- read.csv("./Aux_MonteCarlo/volatilities_larger2_outliers_w_FIGARCH-N_BR.csv")[, -1]^2
figarch_t_out_w <- read.csv("./Aux_MonteCarlo/volatilities_larger2_outliers_w_FIGARCH-T_BR.csv")[, -1]^2
dcs_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_DCS-N_BR.csv")[, -1]^2
dcs_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_DCS-T_BR.csv")[, -1]^2
sv_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_SV-N_BR.csv")[, -1]^2
sv_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_SV-T_BR.csv")[, -1]^2
ms_n_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_MS-N_BR.csv")[, -1]^2
ms_t_out_w <- read.csv("./MonteCarlo/volatilities_outliers_w_MS-T_BR.csv")[, -1]^2


# Tables Bias
bias_dgp_garch <- cbind(table_bias(garch_n), table_bias(garch_n_out_m), table_bias(garch_n_out_w), table_bias(garch_t), table_bias(garch_t_out_m), table_bias(garch_t_out_w))
bias_dgp_figarch <- cbind(table_bias(figarch_n), table_bias(figarch_n_out_m), table_bias(figarch_n_out_w), table_bias(figarch_t), table_bias(figarch_t_out_m), table_bias(figarch_t_out_w))
bias_dgp_dcs <- cbind(table_bias(dcs_n), table_bias(dcs_n_out_m), table_bias(dcs_n_out_w), table_bias(dcs_t), table_bias(dcs_t_out_m), table_bias(dcs_t_out_w))
bias_dgp_sv <- cbind(table_bias(sv_n), table_bias(sv_n_out_m), table_bias(sv_n_out_w), table_bias(sv_t), table_bias(sv_t_out_m), table_bias(sv_t_out_w))
bias_dgp_ms <- cbind(table_bias(ms_n), table_bias(ms_n_out_m), table_bias(ms_n_out_w), table_bias(ms_t), table_bias(ms_t_out_m), table_bias(ms_t_out_w))

tab <- rbind(bias_dgp_garch, bias_dgp_figarch, bias_dgp_dcs, bias_dgp_sv, bias_dgp_ms)
colnames(tab) <- rep(c("500", "1000", "2500"), 6)
row.names(tab) <- rep(c("GARCH-N", "GARCH-T", "FIGARCH-N", "FIGARCH-T", "DCS-N", "DCS-T", "SV-N", "SV-T", "MS-N", "MS-T"), 5)

kbl(tab, format = "latex", booktabs = TRUE, digits = 4) |>
  add_header_above(c(" " = 1,
    "Uncontaminated series" = 3,
    "One additive outlier one month before" = 3,
    "One additive outlier one week before " = 3,
    "Uncontaminated series" = 3,
    "One additive outlier one month before" = 3,
    "One additive outlier one week before " = 3)) |>
  add_header_above(c(" " = 1, "Normal" = 9, "Student-t" = 9)) |>
  save_kable(file =  "tex_tables/bias_MC_BR.tex", keep_tex = TRUE)
  

# Tables MSE, QLIKE, MAE
set.seed(1234)
for (l in c("MSE", "QLIKE")) {
  for (i in c(0, 0.10, 0.25, 0.5)) {
    table_latex(loss_func = l, w = i) |> save_kable(file = paste0("tex_tables/", l, "_", i, "_BR_75_larger2.tex"), keep_tex = TRUE)
  }
}



