## ---------------------------
##
## Snapshot analyzers  
##
## These functions are adapted from https://github.com/bu-rcs/BU-COVID/tree/master/python/bu_covid/bu_covid
##
## ---------------------------


__all__=['get_snapshots','infection_count_to_df','diagnosed_count_to_df','symptomatic_count_to_df','safe_nan_compare']

import covasim as cv
import sciris as sc
import numpy as np
import pandas as pd
from collections import deque
import copy
import covasim.utils as cvu
import itertools as it 
import numba
import operator 

 

def safe_nan_compare(x,y,op):
    ''' do an operation "op" between x & y. np.nan values
        always return a False.'''
    result = np.full_like(x,False,dtype=np.bool)
    if isinstance(x,np.ndarray) and isinstance(y,np.ndarray):
        good_ind= np.isfinite(x) & np.isfinite(y)
        result[good_ind] = op(x[good_ind],y[good_ind])
    elif isinstance(x,np.ndarray):
        good_ind= np.isfinite(x)
        result[good_ind] = op(x[good_ind],y)
    else:
        good_ind= np.isfinite(y)
        result[good_ind] = op(x,y[good_ind])
    return result 

class Daily_infection_count(cv.Analyzer):
    ''' Analyzer that takes a "snapshot" of infected people. '''
       
    def __init__(self, days, *args, **kwargs):
        super().__init__(**kwargs) # Initialize the Analyzer object
        days = sc.promotetolist(days) # Combine multiple days
        days.extend(args) # Include additional arguments, if present
        self.days      = days # Converted to integer representations
        self.dates     = None # String representations
        self.start_day = None # Store the start date of the simulation
        self.snapshots = None # Store the actual snapshots
        return


    def initialize(self, sim):
        self.start_day = sim['start_day'] # Store the simulation start day
        self.days = cv.interventions.process_days(sim, self.days) # Ensure days are in the right format
        self.dates = [sim.date(day) for day in self.days] # Store as date strings
        self.initialized = True
        self.snapshots = {} # Store the snapshots in a dictionary.
        self.yesterday_quarantine = sim.people.quarantined.copy()
        return


    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            today_diag = np.where(ppl.date_exposed.astype(np.int32) == sim.t)
            if len(today_diag[0]) > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['uid'] = ppl.uid[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['GI'] = np.zeros(len(today_diag[0]),dtype=np.uint8)
                self.snapshots[date]['source'] = np.zeros(len(today_diag[0]),dtype=np.uint32)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag[0]]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['GI'][ind] = sim.t - [item['date'] for item in ppl.infection_log if item['target'] == val][0]
                        self.snapshots[date]['exogenous'][ind] = 0 
                        self.snapshots[date]['source'][ind] = ppl.uid[val]

    def get(self, key=None):
        ''' Retrieve a snapshot from the given key (int, str, or date) '''
        if key is None:
            key = self.days[0]
        day  = cv.misc.day(key, start_day=self.start_day)
        date = cv.misc.date(day, start_date=self.start_day, as_date=False)
        if date in self.snapshots:
            snapshot = self.snapshots[date]
        else:
            dates = ', '.join(list(self.snapshots.keys()))
            errormsg = f'Could not find snapshot date {date} (day {day}): choices are {dates}'
            raise sc.KeyNotFoundError(errormsg)
        return snapshot

 

class Daily_diagnosed_count(Daily_infection_count):
    ''' Snapshot the demographics of anyone who is 
        diagnosed on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            today_diag = cv.true(ppl.date_diagnosed.astype(np.int32) == sim.t)
            if  today_diag.size > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['uid'] = ppl.uid[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(today_diag.size,dtype=np.uint8)
                self.snapshots[date]['source'] = np.zeros(today_diag.size,dtype=np.uint32)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['exogenous'][ind] = 0 
                        self.snapshots[date]['source'][ind] = ppl.uid[val]

class Daily_symptomatic_count(Daily_infection_count):
    ''' Snapshot the demographics of anyone who is 
        symptomatic on any given day '''
    def apply(self,sim):
        ppl = sim.people
        for ind in cv.interventions.find_day(self.days, sim.t):
            date = self.dates[ind]
            today_diag = cv.true(ppl.date_symptomatic.astype(np.int32) == sim.t)
            if  today_diag.size > 0:
                self.snapshots[date] = {}
                self.snapshots[date]['uid'] = ppl.uid[today_diag]
                self.snapshots[date]['exogenous'] = np.ones(today_diag.size,dtype=np.uint8)
                self.snapshots[date]['SI'] = np.zeros(today_diag.size,dtype=np.int8)
                self.snapshots[date]['source'] = np.zeros(today_diag.size,dtype=np.uint32)
                source = [item['source'] for item in ppl.infection_log if item['target'] in today_diag]
                for ind, val in enumerate(source):
                    if val is not None:
                        self.snapshots[date]['exogenous'][ind] = 0 
                        self.snapshots[date]['source'][ind] = ppl.uid[val]
                        if not np.isnan(ppl.date_symptomatic[val]): 
                            self.snapshots[date]['SI'][ind] = sim.t - ppl.date_symptomatic[val]


#%%
def infection_count_to_df(sims_complete):
    ''' The infection count is the index 0 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'uid':[], 'exogenous':[],'GI':[],'source':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][0]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            if j+1<len(BU_infect.dates):
                days = BU_infect.days[j]
                date_adj = BU_infect.dates[j+1]
            # Everything that happened on this day gets saved.
                if date_adj in BU_infect.snapshots:
                    for k in range(BU_infect.snapshots[date_adj]['uid'].shape[0]):
                        data['dates'].append(date)
                        data['days'].append(days)
                        data['uid'].append(BU_infect.snapshots[date_adj]['uid'][k])
                        data['exogenous'].append(BU_infect.snapshots[date_adj]['exogenous'][k])
                        data['GI'].append(BU_infect.snapshots[date_adj]['GI'][k])
                        data['source'].append(BU_infect.snapshots[date_adj]['source'][k])
                        count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)

 
def diagnosed_count_to_df(sims_complete):
    ''' The diagnosed count is the index 1 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'uid':[], 'exogenous':[],'source':[]}    
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][1]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            if j+1<len(BU_infect.dates):
                days = BU_infect.days[j]
                date_adj = BU_infect.dates[j+1]
            # Everything that happened on this day gets saved.
                if date_adj in BU_infect.snapshots:
                    for k in range(BU_infect.snapshots[date_adj]['uid'].shape[0]):
                        data['dates'].append(date)
                        data['days'].append(days)
                        data['uid'].append(BU_infect.snapshots[date_adj]['uid'][k])
                        data['exogenous'].append(BU_infect.snapshots[date_adj]['exogenous'][k])
                        data['source'].append(BU_infect.snapshots[date_adj]['source'][k])
                        count += 1            
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
    
def symptomatic_count_to_df(sims_complete):
    ''' The symptomatic count is the index 2 analyzer.  Convert it to 
        a pandas dataframe '''
    data={'sim_num':[], 'dates':[], 'days':[], 'uid':[], 'exogenous':[],'SI':[],'source':[]}
    for i, sim in enumerate(sims_complete):    
        # Get the snapshot
        BU_infect = sim['analyzers'][2]
        # Loop through the dates and add everything found to data.
        # Not the most elegant code, but it gets the job done.
        count = 0 
        for j,date in enumerate(BU_infect.dates):
            if j+1<len(BU_infect.dates):
                days = BU_infect.days[j]
                date_adj = BU_infect.dates[j+1]
            # Everything that happened on this day gets saved.
                if date_adj in BU_infect.snapshots:
                    for k in range(BU_infect.snapshots[date_adj]['uid'].shape[0]):
                        data['dates'].append(date)
                        data['days'].append(days)
                        data['uid'].append(BU_infect.snapshots[date_adj]['uid'][k])
                        data['exogenous'].append(BU_infect.snapshots[date_adj]['exogenous'][k])
                        data['SI'].append(BU_infect.snapshots[date_adj]['SI'][k])
                        data['source'].append(BU_infect.snapshots[date_adj]['source'][k])
                        count += 1
        data['sim_num'] += count * [i]
    return pd.DataFrame(data=data)
 
    
def get_snapshots(num_days):
    ''' Return a list of snapshots to be used with the simulations.  The order here
        is specific and must match that in snapshots_to_df '''
    day_lst = list(range(num_days))
    return [Daily_infection_count(day_lst),
            Daily_diagnosed_count(day_lst),
            Daily_symptomatic_count(day_lst)]
    
