import unittest
import behaviour as bh
import pathosim as inf 
import numpy as np

class test_contact(unittest.TestCase):

    def test_firstDetection(self):

        pars = dict(
            pop_size=1000,
            pop_type='behaviour_module',
            rand_seed = 0,
            start_day = '2020-01-03',
            end_day = '2020-02-13',
            enable_surveillance = True,
            enable_random_testing = True,
            surveillance_viral_threshold = 5,
            surveillance_percentile_threshold = 1,
            enable_syndromic_testing = True,
        )  

        MyNewPathogen = inf.SARS_COV_2(500)
        sim = inf.Sim(pars=pars, pathogens = [MyNewPathogen], verbose = 0)
        simtest = sim.run_as_generator()
        
        for t in simtest:
            if len(sim.viral_loads) == 0:
                return False
            count = sum(1 for load in sim.viral_loads if load > sim.viral_load_threshold)
            percentage = count / len(sim.viral_loads) * 100
            expected_result = percentage >= sim.percent_threshold
            actual_result = sim.percentage >= sim.percent_threshold
            self.assertEqual(expected_result, actual_result)
            
if __name__ == '__main__':
    unittest.main()
       