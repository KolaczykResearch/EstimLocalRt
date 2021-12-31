# EstimLocalRt
This repository contains code supporting the text "Estimation of local time-varying reproduction numbers in noisy surveillance data"



## R

* Create csv file for daily infection counts in Hong Kong: data_hongkong.R
* Create csv file for daily infection counts in Victoria, Australia: data_australia.R
* Bayesian hierarchical model: model_bayesian.R
* Miscellaneous functions used for this paper: Rt_misc.R
* Estimation of local time-varying reproduction numbers for simulated epidemics: simulation_Rt_estimation.R
* Estimation of local time-varying reproduction numbers for real epidemics: application_Rt_estimation.R
* Figure 2: simulation_diagnosed_counts.R
* Figure 3 and Figure 4: simulation_Rt_plot.R
* Figure 5: application_plot.R

## Python

* Library that outputs results from Covasim: Covasim_Rt
* Simulate epidemic in Hong Kong: simulation_epidemic_hongkong.py
* Simulate epidemic in Victoria, Australia: simulation_epidemic_australia.py

## Data

### COVID-19 data in Hong Kong
Adam, Dillon, et al. "Clustering and superspreading potential of SARS-CoV-2 infections in Hong Kong." Nat Med 26, 1714–1719 (2020). 

### COVID-19 data in Victoria, Australia
Seemann, Torsten, et al. "Tracking the COVID-19 pandemic in Australia using genomics." Nature communications 11.1 (2020): 1-9.
