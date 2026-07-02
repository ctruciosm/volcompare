################################################################################
#####                  Tables Empirical Application                        #####
################################################################################
library(dplyr)
library(lubridate)
library(modelconf)
library(tidyr)
library(wrMisc)
library(readxl)
library(kableExtra)
source("Utils_GARCH-GAS-SV.R")


n_ins <- 2500
rm <- "rv"

if (n_ins == 2500){
  datas_oos <- read.csv("./Empirical_Application/data_oos_2500.csv")[, -1]
  garch_n <- read.csv("./Empirical_Application/garch_n_fore_s_2500.csv")[, -1]^2
  garch_t <- read.csv("./Empirical_Application/garch_t_fore_s_2500.csv")[, -1]^2
  figarch_n <- read.csv("./Empirical_Application/figarch_n_fore_s_2500.csv")[, -1]^2
  figarch_t <- read.csv("./Empirical_Application/figarch_t_fore_s_2500.csv")[, -1]^2
  gas_n <- read.csv("./Empirical_Application/gas_n_fore_s_2500.csv")[, -1]^2
  gas_t <- read.csv("./Empirical_Application/gas_t_fore_s_2500.csv")[, -1]^2
  ms_n <- read.csv("./Empirical_Application/ms_n_fore_s_2500.csv")[, -1]
  ms_t <- read.csv("./Empirical_Application/ms_t_fore_s_2500.csv")[, -1]^2
  sv_n <- read.csv("./Empirical_Application/sv_n_fore_s_2500.csv")[, -1]^2
  sv_t <- read.csv("./Empirical_Application/sv_t_fore_s_2500.csv")[, -1]^2
  if (rm == "rv") {
    rv_oos <- read_xlsx("./Data/capire_realised_measures_5min.xlsx", sheet = "RV_5", col_types = c("text", rep("numeric", 30)), na = c("-", " ", "NA")) |> 
      mutate(Date = dmy(Date)) |> filter(Date >= "2019-12-09")
    #rv_oos <- read.csv("./Data/volare_realized_variance_stocks.csv", header = TRUE) |> 
    #  mutate(Date = ymd(date)) |> 
    #  mutate(rv5 = 10000 * rv5) |> 
    #  select(Date, symbol, rv5) |> filter(Date >= "2019-12-09", Date <= "2024-12-31") |> 
    #  pivot_wider(names_from = symbol, values_from = rv5)
  } else {
    rv_oos <- read_xlsx("./Data/capire_realised_measures_5min.xlsx", sheet = "BPV_5", col_types = c("text", rep("numeric", 30)), na = c("-", " ", "NA"))  |> 
      mutate(Date = dmy(Date)) |> filter(Date >= "2019-12-09")
  }
}


m_mse <- matrix(0, ncol = 10, nrow = ncol(garch_n))
m_qlike <- matrix(0, ncol = 10, nrow = ncol(garch_n))

for (i in 1:ncol(garch_n)) {
  m_mse[i, ] <- c(mean(loss_mse(rv_oos[, colnames(garch_n)[i]], garch_n[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], garch_t[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], figarch_n[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], figarch_t[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], gas_n[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], gas_t[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], ms_n[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], ms_t[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], sv_n[, colnames(garch_n)[i]])[[1]]),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], sv_t[, colnames(garch_n)[i]])[[1]]))
  m_qlike[i, ] <- c(mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], garch_n[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], garch_t[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], figarch_n[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], figarch_t[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], gas_n[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], gas_t[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], ms_n[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], ms_t[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], sv_n[, colnames(garch_n)[i]])[[1]]),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], sv_t[, colnames(garch_n)[i]])[[1]]))
}


col_names <- c("GARCH-N", "GARCH-T", "FIGARCH-N", "FIGARCH-T", "GAS-N", "GAS-T", "MS-N", "MS-T", "SV-N", "SV-T")
colnames(m_mse) <- col_names
colnames(m_qlike) <- col_names
row.names(m_mse) <- colnames(garch_n)
row.names(m_qlike) <- colnames(garch_n)

################################################################################
#####                              MSC                                     #####
################################################################################

stock_names <- colnames(garch_n)
mcs_mse <- mcs_qlike <- matrix(0, ncol = 10, nrow = 29)
set.seed(123)
for (i in 1:ncol(garch_n)) {
  aux_mse <- cbind(loss_mse(rv_oos[, stock_names[i]], garch_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], garch_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], figarch_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], figarch_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], gas_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], gas_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], ms_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], ms_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], sv_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], sv_t[, stock_names[i]]))
  
  aux_qlike <- cbind(loss_qlike(rv_oos[, stock_names[i]], garch_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], garch_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], figarch_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], figarch_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], gas_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], gas_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], ms_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], ms_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], sv_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], sv_t[, stock_names[i]]))
  
  colnames(aux_mse) <- col_names
  colnames(aux_qlike) <- col_names
  
  mcs_mse[i, estMCS.quick(aux_mse, test = "t.range", B = 10000, l = 21, alpha = 0.25)] <- 1
  mcs_qlike[i, estMCS.quick(aux_qlike, test = "t.range", B = 10000, l = 21, alpha = 0.25)] <- 1
}

################################################################################
#####                            Latex                                     #####
################################################################################

table_metrics <- m_mse
table_mcs <- mcs_mse
table_latex <- table_metrics|> kbl(format = "latex", booktabs = TRUE, digits = 3)
  
for (j in 1:dim(table_metrics)[2]) {
  bg <- ifelse(table_mcs[, j] == 1, "gray!45", "white")
  table_latex <- table_latex |> column_spec(j + 1, background = bg)
}
table_latex |> save_kable(file = paste0("Emp_App_", n_ins, "_", rm, "mse", ".tex"), keep_tex = TRUE)


table_metrics <- m_qlike
table_mcs <- mcs_qlike
table_latex <- table_metrics|> kbl(format = "latex", booktabs = TRUE, digits = 3)

for (j in 1:dim(table_metrics)[2]) {
  bg <- ifelse(table_mcs[, j] == 1, "gray!45", "white")
  table_latex <- table_latex |> column_spec(j + 1, background = bg)
}
table_latex |> save_kable(file = paste0("Emp_App_", n_ins, "_", rm, "qlike", ".tex"), keep_tex = TRUE)







