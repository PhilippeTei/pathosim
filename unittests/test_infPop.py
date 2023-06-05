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
            'n_days': 10,
            'pop_infected': 50,
            'rand_seed' : 0
        }

class test_infPop(unittest.TestCase):

    def test_different_pop_types(self):

        baseline_contacts = sc.loadobj(f'{full_path}/test_infPop_contacts.baseline') 
        baseline_pars = sc.loadobj(f'{full_path}/test_infPop_pars.baseline') 

        if(baseline_contacts == None):
            raise FileNotFoundError 
        if(baseline_pars == None):
            raise FileNotFoundError 

        poptype = 'behaviour_module'

        pop = sc.loadobj(f'{full_path}/test_infPop_bhPopulation.ppl')
        sim = inf.Sim(simPars, pop_type = poptype, people = pop, verbose = 0)
        sim.init_people()  
           
        #check all keys are present
        for k in baseline_pars:
            self.assertEqual(True if k in sim.people.pars else False, True)

        #check some parameters
        self.assertEqual(True if np.array_equal(sim.people.pars['prognoses']['symp_probs'], baseline_pars['prognoses']['symp_probs']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people.pars['prognoses']['severe_probs'], baseline_pars['prognoses']['severe_probs']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people.pars['prognoses']['crit_probs'], baseline_pars['prognoses']['crit_probs']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people.pars['prognoses']['death_probs'], baseline_pars['prognoses']['death_probs']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people.pars['prognoses']['comorbidities'], baseline_pars['prognoses']['comorbidities']) else False, True)

        self.assertEqual(True if np.array_equal(sim.people['contacts']['h']['p1'], baseline_contacts['h']['p1']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['h']['p2'], baseline_contacts['h']['p2']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['s']['p1'], baseline_contacts['s']['p1']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['s']['p2'], baseline_contacts['s']['p2']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['w']['p1'], baseline_contacts['w']['p1']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['w']['p2'], baseline_contacts['w']['p2']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['c']['p1'], baseline_contacts['c']['p1']) else False, True)
        self.assertEqual(True if np.array_equal(sim.people['contacts']['c']['p2'], baseline_contacts['c']['p2']) else False, True)
                 

def generate_baseline():
    import behaviour as bh
    pop_size = 1000
    seed = 0 
    bh_pars = sc.objdict(n = pop_size,
                            rand_seed = seed,
                            country_location = "england")
        
    pop = bh.BehaviourModel(bh_pars)
    
    pop.save(f'{full_path}/test_infPop_bhPopulation.ppl')
     
    sim = inf.Sim(simPars,pop_type ='behaviour_module', people = pop, verbose = 0)
    sim.init_people() 
    sc.saveobj(f'{full_path}/test_infPop_contacts.baseline', sim.people.contacts)
    sc.saveobj(f'{full_path}/test_infPop_pars.baseline', sim.people.pars)

    

if __name__ == '__main__': 
    if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
    else:
       unittest.main()