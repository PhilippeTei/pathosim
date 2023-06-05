
import unittest
import sciris as sc
import behaviour as bh

import os

gen_baseline = False 

absolute_path = os.path.dirname(__file__)
relative_path = "baselines"
full_path = os.path.join(absolute_path, relative_path)

class test_bhPop(unittest.TestCase):
    

    def test_simplePop(self):
        pop_size = 20000
        seed = 0

        bh_pars = sc.objdict(n = pop_size,
                                rand_seed = seed,
                                country_location = "england")
        
        pop = bh.BehaviourModel(bh_pars)

        baseline = sc.loadobj(f'{full_path}/test_bhPop.baseline')

        if(baseline == None):
            raise FileNotFoundError

        self.assertEqual(pop.popdict, baseline)
         
 
        

def generate_baseline():
    pop_size = 20000
    seed = 0

    bh_pars = sc.objdict(n = pop_size,
                            rand_seed = seed,
                            country_location = "england")
        
    pop = bh.BehaviourModel(bh_pars)
    
    sc.save(f'{full_path}/test_bhPop.baseline',pop.popdict)

if __name__ == '__main__':

   if gen_baseline:
        print("Generating baseline file")
        generate_baseline()
   else:
       unittest.main()
