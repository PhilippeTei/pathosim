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
vx = inf.vaccinate_num('pfizer',200, False,None, None)
vxbst = inf.vaccinate_num('pfizer',150, True,None, None)
 
cb = inf.change_beta(days=40, changes=0.5) 

 
nc = inf.notify_contacts(trace_probs={'h': 1, 's': 0.2, 'w': 0.2, 'c': 0.03},
                                  trace_time={'h': 0, 's': 1, 'w': 1, 'c': 2})

RATcapacity = np.ones(50+1)*int(0.1*1000)
PCRcapacity = np.ones(50+1)*int(0.01*1000)

RAT_di_criteria = ['cont_conf', 'cont_vuln', 'cont_sx', 'sx', 'other']
RAT_sc_criteria = ['work']
PCR_di_criteria = ['sx', 'other']
PCR_sc_criteria = ['work']

RATpars = dict(capacity=RATcapacity, di_criteria=RAT_di_criteria, sc_criteria=RAT_sc_criteria, test_mode='LOD')
PCRpars = dict(capacity=PCRcapacity, di_criteria=PCR_di_criteria, sc_criteria=PCR_sc_criteria, test_mode='LOD')

test_objs = dict(PCR_disc=PCRpars, RAT_disc=RATpars)


pars1 = dict(
        use_waning    = True,          
        enable_vl = True,              
        pop_size      = 5000,          
        pop_infected  = 200,           
        pop_type      = 'random',      
        n_days        = 100,            
        verbose       = 0,             
        rand_seed     = 0,             
        interventions = [vx, vxbst],           
    )                                  
pars2 = dict(
        use_waning    = True,          
        enable_vl = True,              
        pop_size      = 5000,          
        pop_infected  = 50,           
        pop_type      = 'random',      
        n_days        = 100,            
        verbose       = 0,             
        rand_seed     = 1,             
        interventions = [cb],           
    )        
pars3 = dict(
        use_waning    = True,          
        enable_vl = True,              
        pop_size      = 1000,          
        pop_infected  = 50,           
        pop_type      = 'behaviour_module',      
        n_days        = 50,            
        verbose       = 0,             
        rand_seed     = 3,             
        interventions = [nc],   
        enable_testobjs = True,
        testobjs = test_objs
    )        
   
result_keys = [ 
    'n_infectious',  
    'n_symptomatic', 
    'n_severe',       
    'n_recovered',   
    'n_dead']


class test_covSimWithInterventions(unittest.TestCase):
    
    def test_vaccinateNumAndBoosters(self):
        sim = inf.Sim(pars1)
        sim.run()
        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_vaccinateNumAndBoosters_{k}.baseline')  
            self.assertEqual(True if expected_result==sim.results[k] else False, True)

    def test_changeBeta(self):
        sim = inf.Sim(pars2)
        sim.run()
        for k in result_keys:
            expected_result = sc.loadobj(f'{full_path}/test_changeBeta_{k}.baseline')  
            self.assertEqual(True if expected_result==sim.results[k] else False, True)
    '''
    def test_contactTracing(self):
        pop = sc.loadobj(f'{full_path}/test_infPop_bhPopulation.ppl')
        sim = inf.Sim(pars3, people = pop)
        sim.run()
        
        expected_known_contact = sc.loadobj(f'{full_path}/test_contactTracing_kc.baseline')   
        expected_date_known_contact = sc.loadobj(f'{full_path}/test_contactTracing_dkc.baseline')  
         
        for i in range(len(expected_known_contact)):
            self.assertEqual(expected_known_contact[i], sim.people.known_contact[i])'''
                 

def generate_baseline():
     
    sim = inf.Sim(pars1)
    sim.run()
    for k in result_keys:
        sc.saveobj(f'{full_path}/test_vaccinateNumAndBoosters_{k}.baseline', sim.results[k]) 

    sim = inf.Sim(pars2)
    sim.run()
    for k in result_keys:
        sc.saveobj(f'{full_path}/test_changeBeta_{k}.baseline', sim.results[k]) 

    '''
    pop = sc.loadobj(f'{full_path}/test_infPop_bhPopulation.ppl')
    sim = inf.Sim(pars3, people = pop)
    sim.run()
          
    sc.saveobj(f'{full_path}/test_contactTracing_kc.baseline',sim.people.known_contact)   
    sc.saveobj(f'{full_path}/test_contactTracing_dkc.baseline',sim.people.date_known_contact)   '''

    return

    

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
        
       unittest.main()