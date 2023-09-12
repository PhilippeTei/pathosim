 
import pathosim as inf
# Reload if it's changed. 
import importlib
importlib.reload(inf)
import numpy as np
from behaviour import RegionalBehaviourModel
importlib.reload(inf)
import matplotlib.pyplot as plt
from behaviour import BehaviourModel
import sciris as sc
import os

""" Parameters for Behaviour Module """
## Specify demographic parameters. Given as a list. ##
reg_sizes = [4000] * 5
reg_names = ['A', 'B', 'C', 'D', 'E']
reg_names = [f'Region {name}' for name in reg_names]

reg_params = []
for reg in reg_names:
    cur_params = dict(
        name=reg,
        n=reg_sizes[reg_names.index(reg)],
        com_contacts=10,
        country_location='england'
    )
    reg_params.append(cur_params)

## Specify work mixing parameters. 20% leave, to other places uniformly. ##

params_work_mixing = {}

for reg in reg_names:
    # People leave uniformly to other regions.
    other_regions = [x for x in reg_names if x != reg]
    dests = {}

    for other_reg in other_regions:
        dests[other_reg] = 1 / len(other_regions)
    
    # Ensure that dests sums to 1. 
    dests[other_regions[-1]] = 1 - sum([dests[x] for x in other_regions[:-1]])

    params_work_mixing[reg] = {"leaving": 0.2, "dests": dests}

"""
Aside: a simple hard-coded 3 region example looks like:

dests_a = {'city_b':0.6, 'city_c':0.4}
dests_b = {'city_a':0.9, 'city_c':0.1}
dests_c = {'city_a':0.9, 'city_b': 0.1}

params_work_mixing = dict(city_a = {"leaving":0.05, "dests":dests_a},
                            city_b =    {"leaving":0.4, "dests":dests_b},
                            city_c =  {"leaving":0.3, "dests":dests_c})

"""

## Specify community mixing parameters. 20% community-based contacts are 
# with people from other communities, uniformly. ##
with_others = 0.2

params_com_mixing = (with_others/(len(reg_names)-1))*np.ones((len(reg_names), len(reg_names)))

# Populate the diagonal with 1 - with_others.
for i in range(params_com_mixing.shape[0]):
    params_com_mixing[i][i] = 1-with_others

# Create pop object
pop = RegionalBehaviourModel(com_mixing=True, work_mixing=True, params_com_mixing=params_com_mixing, 
                params_work_mixing=params_work_mixing, all_reg_params=reg_params)

## Specify demographic parameters. Given as a list. ## 
region_seed_infections = dict(
    {'Region A': 20, 'Region B': 15, 'Region C': 10, 'Region D': 5, 'Region E': 20}
)

multi_reg_params = dict( 
    reg_params=pop.reg_pars # Passes various parameters for behaviour module to Covasim. 
)

pars = dict(
    pop_size=20000, 
    pop_type='behaviour_module',
    n_days = 200,
    rand_seed = 0
)
cov = inf.SARS_COV_2(region_seed_infections)
pars_mr = pars 
pars_mr['multiregion'] = multi_reg_params
pars_mr['enable_multiregion'] = True
sim_mr = inf.Sim(pars=pars_mr, people=pop, pathogens = cov)
sim_mr.run()  


 