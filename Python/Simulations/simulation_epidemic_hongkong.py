#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import covasim as cv
import sciris as sc
import os 
import numpy as np


# This ***MUST*** be imported after covasim.
import Covasim_Rt as rt

# Name of this simulation for output files 
sim_name = 'hongkong'

# =============================================================================
 
# File for daily imported counts 
pop_imported_path = '../../../Data/HongKong/imported_counts_adjusted.csv' 

# Set a destination for output files 
plot_dir = '../../../Results/Simulation_Epidemics/'
os.makedirs(plot_dir, exist_ok = True) 

# Set the number of parallel runs 
n_runs = 1000

# =============================================================================
#  Build intervention
# =============================================================================
# The simulation starts
start_day = '2020-01-18'
end_day = '2020-04-26'

# Total number of simulation days
num_days = (rt.make_dt(end_day) - rt.make_dt(start_day)).days +1
pop_size = 1e5
beta_val = 0.016
pars = dict(pop_size = pop_size,
            pop_infected = 0, # Initial number of infected people
            beta = beta_val,      # Base transmission rate
            start_day = start_day,
            end_day = end_day,
            quar_factor = 0,
            iso_factor = 0,
            asymp_factor=0.5,
            n_imports=0,
            rescale=False,
            contacts = 20) 
            
import_exogenous = np.genfromtxt(pop_imported_path, delimiter=',')  

base_interventions = []
daily_n_imports = np.genfromtxt(pop_imported_path, delimiter=',')   
n_imports_days = np.arange(num_days)
# Details of these functions are available at https://github.com/InstituteforDiseaseModeling/covasim/tree/master/covasim
base_interventions.append( cv.dynamic_pars({'n_imports':{'days':n_imports_days, 'vals':daily_n_imports}}))
base_interventions.append(cv.test_num(daily_tests=round(pop_size/3), start_day = 0, test_delay = 1, quar_test=1, symp_test = 1,sensitivity= 0.9))
base_interventions.append(cv.clip_edges([60, 70], [0.8, 1], layers='a'))
base_interventions.append(cv.change_beta([60, 70], [0.8, 1], layers='a'))
base_interventions.append(cv.clip_edges([70, 100], [0.4, 1], layers='a'))
base_interventions.append(cv.change_beta([70, 100], [0.4, 1], layers='a'))  
interventions = base_interventions.copy()

# =============================================================================
#  Run the simulations
# =============================================================================

verbose = False # keep the simulations quiet
if n_runs == 1:
    verbose = True  #unles there's just 1.

#%% Create the list of daily snapshot objects
analyzers = rt.get_snapshots(num_days)

#%% Run the simulations
sim=cv.Sim(pars=pars, verbose = verbose, interventions = interventions, analyzers=analyzers)
sims_complete = rt.parallel_run_sims(sim, n_runs = n_runs, n_cores = rt.get_n_cores())
msim = cv.MultiSim(sims_complete)
msim.reduce()


print("Creating dataframes") 
#%% Get the snapshots dataframe
infection_df=rt.infection_count_to_df(sims_complete)
diagnosed_df = rt.diagnosed_count_to_df(sims_complete)
symptomatic_df = rt.symptomatic_count_to_df(sims_complete)
sim_results_df = rt.sim_results_to_df(sims_complete)


# =============================================================================
#  Output
# =============================================================================

with  open(os.path.join(plot_dir,'sim_results_%s.csv' % sim_name), 'w') as f:
    sim_results_df.to_csv(f)

with  open(os.path.join(plot_dir,'infection_counts_%s.csv' % sim_name), 'w') as f:
    infection_df.to_csv(f)

with  open(os.path.join(plot_dir,'diagnosed_counts_%s.csv' % sim_name), 'w') as f:
    diagnosed_df.to_csv(f)      

with  open(os.path.join(plot_dir,'symptomatic_counts_%s.csv' % sim_name), 'w') as f:
    symptomatic_df.to_csv(f)     

 