import behaviour as bh
import pathosim as inf 
import numpy as np
import pandas as pd
import os

pars = dict(
pop_size=5000,
pop_type='random',
n_days = 100,
rand_seed = 0,
start_day = '2020-01-30',
enable_surveillance = True,
enable_syndromic_testing = True,
)

Coronavirus = inf.SARS_COV_2(5)
sim = inf.Sim(pars=pars, pathogens = [Coronavirus], verbose = 1)
sim.run()
print(sim.total_runs)
print(sim.total_costs)

#for i in range(3):
    #Coronavirus = inf.SARS_COV_2(100)
    #sim = inf.Sim(pars=pars, pathogens = [Coronavirus], verbose = 1)
    #sim.run()

# simtest = sim.run_as_generator()

# for t in simtest:
#     print(sim.detection)