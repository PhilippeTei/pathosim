import covasim as cv
import pandas as pd
import sciris as sc
import os
import multiprocessing

ukdata = pd.read_csv(r"C:\Simdata\ukearlydata.csv")

pars = sc.objdict(
    pop_size = 100000,
    scaled_pop = 3503238.747745823,
    pop_infected = 182,
    pop_type = 'hybrid',
    location = "United Kingdom of Great Britain and Northern Ireland", 
    start_day = '2020-01-30',
    end_day   = '2020-09-30',
    beta = 0.031064847611335867,
    verbose = 0,
    interventions = [cv.test_num(daily_tests='data'), cv.change_beta(days=['2020-02-26', '2020-07-26'], changes=[0.22, 0.7])],
    rel_death_prob = 2.8264951851502245,
    rel_symp_prob = 0.825205571316789,
    rel_severe_prob = 2.305903581404782,
    rel_crit_prob = 10.7178410896736656,
    n_imports = 100
)

sim = cv.Sim(pars=pars, datafile=ukdata)
sim.run()
fit = sim.compute_fit()
fit.summarize()
fit.plot()
