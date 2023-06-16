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
            'rand_seed' : 0,
            'verbose':0 
        }

result_keys = [ 
    'n_infectious',  
    'n_symptomatic', 
    'n_severe',       
    'n_recovered',   
    'n_dead']

class test_covVariantsSim(unittest.TestCase):
    
    def test_oneVariant(self): 
        
        pathogen = inf.SARS_COV_2(20)
        pathogen.add_existing_variant('alpha', days=5, n_imports=10) 
         
        sim = inf.Sim(simPars, pathogens = [pathogen]) 
        sim.run()
         
        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_covVariantsSim_oneVariant_{k}.baseline')  
            self.assertEqual(True if expected_result==sim.results[k] else False, True)
            #if expected_result!=sim.results[k]:
              #  print(sim.result[k])
              #  print(expected_result)
              #  print('---------------')
    
    def test_multiVariants(self): 
        pathogen = inf.SARS_COV_2(20)

        pathogen.add_existing_variant('alpha', days=0, n_imports=10)
        pathogen.add_existing_variant('beta', days=10, n_imports=10)
        pathogen.add_existing_variant('omicron', days=30, n_imports=15)
        pathogen.add_custom_variant(label='3x more transmissible', rel_beta = 3.0, days=60, n_imports=20)
         
        sim = inf.Sim(simPars, pathogens = [pathogen]) 
        sim.run()

        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_covVariantsSim_multiVariant_{k}.baseline')  
            if not expected_result==sim.results[k]:
                print(expected_result)
                print("___________")
                print(sim.results[k])
            self.assertEqual(True if expected_result==sim.results[k] else False, True)

                 

def generate_baseline():
    
    pathogen = inf.SARS_COV_2(20)
    
    pathogen.add_existing_variant('alpha', days=5, n_imports=10) 
         
    sim = inf.Sim(simPars, pathogens = [pathogen]) 
    sim.run() 
    
    for k in result_keys:
        sc.saveobj(f'{full_path}/test_covVariantsSim_oneVariant_{k}.baseline', sim.results[k]) 
        
    
    pathogen = inf.SARS_COV_2(20)
     
    pathogen.add_existing_variant('alpha', days=0, n_imports=10)
    pathogen.add_existing_variant('beta', days=10, n_imports=10)
    pathogen.add_existing_variant('omicron', days=30, n_imports=15)
    pathogen.add_custom_variant(label='3x more transmissible', rel_beta = 3.0, days=60, n_imports=20)
        
    sim = inf.Sim(simPars, pathogens = [pathogen]) 
    sim.run() 

    for k in result_keys:
        sc.saveobj(f'{full_path}/test_covVariantsSim_multiVariant_{k}.baseline', sim.results[k]) 

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
       unittest.main()