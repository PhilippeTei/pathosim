import unittest 
import infection as inf 
import sciris as sc
import os
import numpy as np
gen_baseline = False 

absolute_path = os.path.dirname(__file__)
relative_path = "baselines"
full_path = os.path.join(absolute_path, relative_path)

simPars = {
            'pop_size': 1000,
            'n_days': 120,
            'pop_infected': 20,
            'rand_seed' : 0,
            'verbose':0
        }

result_keys = [ 
    'n_infectious',   
    'n_severe',       
    'n_recovered',   
    'n_dead']

class test_covVariantsSim(unittest.TestCase):
    
    def test_oneVariant(self): 
        alpha  = inf.variant('alpha', days=5, n_imports=10)  
        
        sim = inf.Sim(simPars, variants=[alpha]) 
        sim.run()

        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_covVariantsSim_oneVariant_{k}.baseline')  
            self.assertEqual(True if expected_result==sim.results[k] else False, True)
    
    def test_multiVariants(self): 
        alpha  = inf.variant('alpha', days=0, n_imports=10)
        beta   = inf.variant('beta', days=10, n_imports=10)
        omicron   = inf.variant('omicron', days=30, n_imports=15)
        custom = inf.variant(label='3x more transmissible', variant={'rel_beta': 3.0}, days=60, n_imports=20)
        
        sim = inf.Sim(simPars, variants=[alpha, beta, omicron, custom]) 
        sim.run()

        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_covVariantsSim_multiVariant_{k}.baseline')  
            self.assertEqual(True if expected_result==sim.results[k] else False, True)

                 

def generate_baseline():
    
    alpha  = inf.variant('alpha', days=5, n_imports=10)
    
    sim = inf.Sim(simPars, variants=[alpha]) 
    sim.run() 
    
    for k in result_keys:
        sc.saveobj(f'{full_path}/test_covVariantsSim_oneVariant_{k}.baseline', sim.results[k]) 
        

    alpha  = inf.variant('alpha', days=0, n_imports=10)
    beta   = inf.variant('beta', days=10, n_imports=10)
    omicron   = inf.variant('omicron', days=30, n_imports=15)
    custom = inf.variant(label='3x more transmissible', variant={'rel_beta': 3.0}, days=60, n_imports=20)
        
    sim = inf.Sim(simPars, variants=[alpha, beta, omicron, custom]) 
    sim.run() 

    for k in result_keys:
        sc.saveobj(f'{full_path}/test_covVariantsSim_multiVariant_{k}.baseline', sim.results[k]) 

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
       unittest.main()