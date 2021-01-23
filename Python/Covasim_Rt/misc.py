## ---------------------------
## Miscellaneous functions used for file output, etc.
## 
## These functions are adapted from https://github.com/bu-rcs/BU-COVID/tree/master/python/bu_covid/bu_covid
##
## ---------------------------

import numpy as np
import datetime
import os 
import pandas as pd
import networkx as nx
import zipfile
import tqdm


import covasim as cv
 

__all__=['make_dt','sim_results_to_df', ]
 
#%%
 
def make_dt(day):
    ''' Make a Python datetime object. day is YYYY-MM-DD'''
    return datetime.datetime.strptime(day,'%Y-%m-%d') 


#%%
def sim_results_to_df(sims_complete):
    ''' Takes a list of completed sims.  Returns
        a Pandas dataframe of all of their results. The sim_num
        column indexes the simulations. '''
    data={'sim_num':[], 'dates':[], 'days':[]}
    # Take the 1st simulation and add all of its results keys to the dictionary
    sim0 = sims_complete[0]
    keys = list(sim0.results.keys())
    # Remove day num and date keys
    keys.pop(keys.index('t'))
    keys.pop(keys.index('date'))
    # Add the rest.
    for k in keys:
        data[k] = []
    
    for i, sim in enumerate(sims_complete):
        # Convert all quantities to python lists to avoid the
        # overhead of ndarray concatentation
        for k in keys:
            data[k] += list(sim.results[k])
        data['days'] += list(sim.results['t'])
        data['dates'] += list(sim.results['date'])
        data['sim_num'] += len(sim.results['date']) * [i]
        
    return pd.DataFrame(data=data)
 

