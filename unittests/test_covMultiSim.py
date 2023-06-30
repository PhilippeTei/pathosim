import unittest 
import pathosim as inf 
import sciris as sc
import os
import numpy as np

gen_baseline = False 

absolute_path = os.path.dirname(__file__)
relative_path = "baselines"
full_path = os.path.join(absolute_path, relative_path)
 

pars = dict(
        use_waning    = True,            
        pop_size      = 5000,           
        pop_type      = 'random',      
        n_days        = 80,            
        verbose       = 0,             
        rand_seed     = 0                     
    )                                  

result_keys = [ 
    'n_infectious',   
    'n_severe',       
    'n_recovered',   
    'n_dead']


class test_covMultiSim(unittest.TestCase):

    def test_multiSim(self): 
     
        betas = np.linspace(0.010, 0.020, 5) # Sweep beta from 0.01 to 0.02 with 5 values
        sims = []
        for beta in betas:
            sim = inf.Sim(pars, beta=beta, label=f'Beta = {beta}')
            sims.append(sim)
        msim = inf.MultiSim(sims, verbose = 0)
        msim.run()
        msim.reduce()

        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_covMultiSim{k}.baseline')  
            self.assertEqual(True if expected_result==msim.results[0][k] else False, True)


             
def generate_baseline():
     
    betas = np.linspace(0.010, 0.020, 5) # Sweep beta from 0.01 to 0.02 with 5 values
    sims = []
    for beta in betas:
        sim = inf.Sim(pars, beta=beta, label=f'Beta = {beta}')
        sims.append(sim)
    msim = inf.MultiSim(sims, verbose = 0)
    msim.run()
    msim.reduce()

    for k in result_keys:
        sc.saveobj(f'{full_path}/test_covMultiSim{k}.baseline', msim.results[0][k])   

    return

    

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
       unittest.main()