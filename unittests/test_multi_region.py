import unittest 
import pathosim as inf   
import numpy as np
from behaviour import RegionalBehaviourModel 
import matplotlib.pyplot as plt
from behaviour import BehaviourModel
import sciris as sc
import os

gen_baseline = False

absolute_path = os.path.dirname(__file__)
relative_path = "baselines"
full_path = os.path.join(absolute_path, relative_path)

simPars = {
            'pop_size': 1000,
            'n_days': 120, 
            'rand_seed' : 0,
            'verbose':0 
        }

keys_to_check = ['n_infectious', 'n_dead', 'n_symptomatic'] 
regs = ['A','B','C','D','E']

class test_multi_region(unittest.TestCase):
    
    def test_multi_reg(self): 
        
        pop = sc.load(f'{full_path}/test_multi_reg_pop.pop')

        ## Specify demographic parameters. Given as a list. ## 
        region_seed_infections = dict(
            {'Region A': 20, 'Region B': 15, 'Region C': 10, 'Region D': 5, 'Region E': 20}
        )

        multi_reg_params = dict( 
            reg_params=pop.reg_pars # Passes various parameters for behaviour module to Covasim. 
        )

        pars = dict(
            pop_size=25000, 
            pop_type='behaviour_module',
            n_days = 200,
            rand_seed = 0
        )
        cov = inf.SARS_COV_2(region_seed_infections)
        pars_mr = pars 
        pars_mr['multiregion'] = multi_reg_params
        pars_mr['enable_multiregion'] = True
        sim_mr = inf.Sim(pars=pars_mr, people=pop, pathogens = [cov], verbose = 0)
        sim_mr.run() 
 
        for key in keys_to_check:  
            index = 0
            for reg_name in region_seed_infections.keys():
                s = 0
                for i in sim_mr.results[0][f'{reg_name}_{key}']: 
                    s += i  
                exp = sc.load(f'{full_path}/test_multi_region_Region_{regs[index]}_{key}.baseline')   
                index += 1 
                assert s == exp

        for key in keys_to_check:   
            s = sim_mr.results[0][f'{key}']
            exp = sc.loadobj(f'{full_path}/test_multi_region_country_{key}.baseline') 
            assert np.array_equal(s, exp) == True
                 

def generate_baseline():

    pop = sc.load(f'{full_path}/test_multi_reg_pop.pop')

    ## Specify demographic parameters. Given as a list. ## 
    region_seed_infections = dict(
        {'Region A': 20, 'Region B': 15, 'Region C': 10, 'Region D': 5, 'Region E': 20}
    )

    multi_reg_params = dict( 
        reg_params=pop.reg_pars # Passes various parameters for behaviour module to Covasim. 
    )

    pars = dict(
        pop_size=25000, 
        pop_type='behaviour_module',
        n_days = 200,
        rand_seed = 0
    )
    cov = inf.SARS_COV_2(region_seed_infections)
    pars_mr = pars 
    pars_mr['multiregion'] = multi_reg_params
    pars_mr['enable_multiregion'] = True
    sim_mr = inf.Sim(pars=pars_mr, people=pop, pathogens = [cov], verbose = 0)
    sim_mr.run() 
    keys_to_check = ['n_infectious', 'n_dead', 'n_symptomatic']

    for key in keys_to_check:  
        index = 0
        for reg_name in region_seed_infections.keys():
            s = 0
            for i in sim_mr.results[0][f'{reg_name}_{key}']: 
                s += i  
            sc.saveobj(f'{full_path}/test_multi_region_Region_{regs[index]}_{key}.baseline', s)   
            index += 1 
          
    for key in keys_to_check:    
        sc.saveobj(f'{full_path}/test_multi_region_country_{key}.baseline', sim_mr.results[0][f'{key}'])  
 

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
       unittest.main()