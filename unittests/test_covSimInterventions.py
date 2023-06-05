import unittest 
import infection as inf 
import sciris as sc
import os
import numpy as np

gen_baseline = False 

absolute_path = os.path.dirname(__file__)
relative_path = "baselines"
full_path = os.path.join(absolute_path, relative_path)

# Define the interventions
tp = inf.test_prob(start_day=20, symp_prob=0.1, asymp_prob=0.01)
vx = inf.vaccinate_prob('pfizer', days=30, prob=0.1)
cb = inf.change_beta(days=40, changes=0.5)
ct = inf.contact_tracing(trace_probs=0.3, start_day=50)

pars1 = dict(
        use_waning    = True,          
        enable_vl = True,              
        pop_size      = 5000,          
        pop_infected  = 50,           
        pop_type      = 'random',      
        n_days        = 50,            
        verbose       = 0,             
        rand_seed     = 2,             
        interventions = [tp, vx, cb, ct],           
    )                                  

result_keys = [
    'n_susceptible',
    'n_infectious',  
    'n_symptomatic', 
    'n_severe',      
    'n_critical',    
    'n_recovered',   
    'n_dead']


class test_covSimWithInterventions(unittest.TestCase):

    def test_interventions(self): 
        sim = inf.Sim(pars1)
        sim.run()
        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_covSimWithInterventions_{k}.baseline')  
            self.assertEqual(True if expected_result==sim.results[k] else False, True)


                 

def generate_baseline():
     
    sim = inf.Sim(pars1)
    sim.run()
    for k in result_keys:
        sc.saveobj(f'{full_path}/test_covSimWithInterventions_{k}.baseline', sim.results[k]) 

    return

    

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
       unittest.main()