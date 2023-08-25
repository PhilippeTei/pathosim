import behaviour as bh
import pathosim as inf 
import numpy as np
import copy as cp
 
pop_size = 2000

MyNewPathogen = inf.SARS_COV_2(600)   
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 20, verbose = 1, enable_surveillance = True, enable_syndromic_testing = True)
#sim.run()
#sim.plot()

sim.run_as_generator()

for t in sim:
    print(t)