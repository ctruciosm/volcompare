# Replication Codes

This repository contains the codes used to replicate the results in Marques and Trucíos (2026).

### Empirical Application

- `Empirical_Application.R` computes one-step-ahead conditional variances and one-step-ahead VaR and ES for the daily returns series of the constituents of the Dow Jones Industrial Average Index.
- `Tables_App.R` Performs the Model Confidence Set procedure for the empirical application results (also generates the Tables 5 in the main manuscript and 10 in the Supplementary Material).
- `Tables_App_VaR_ES.R` Performs the calibration tests and apply the MCS to the scoringh functions for the VaR and ES (also generates Table 6)

> Five-minute realized variances are freely available from the [CaPiRe](https://capire.stat.unipd.it/) database. Daily returns were obtained from Economatica.

### Monte Carlo Simulation

- `MonteCarloSimulation.R` runs the one-step-ahead forecasting experiment. To use the code, modify the parameters accordingly, or execute it in batch mode using, for instance, the following command:  
  `R CMD BATCH "--args GARCH-N BR" MonteCarloSimulations.R MonteCarlo_GARCH-N_BR.txt &`  
  (You can change `BR` to `US`, `FALSE` to `TRUE`.)
- `Tables_MonteCarlo.R` reproduces the results shown in Tables 2 and 3 of the main manuscript as well as Tables 4 - 9 in the Supplementary Material
- `Aux_MonteCarlo.R` reproduces Table 1 in the Supplementary Material.
- `MonteCarloSimulations_Larger_Sample` reproduces Table 2 in the Supplementary Material.

### Auxiliary Functions

- `DGPs.R` defines the data-generating processes used in the simulations.
- `Utils_GARCH-GAS-SV.R` and `utils.cpp` contain additional functions for model estimation and forecasting.
- `Descriptive_Statistics` displays the descriptive statistics in Table 4.


## References

Marques, F. and Trucíos C. (2026). *"Daily Volatility Forecasting with Off-the-Shelf Models: A Comparative Study Under Stress".* Submitted
