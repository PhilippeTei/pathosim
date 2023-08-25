import unittest
import behaviour as bh
import pathosim as inf 
import numpy as np

class TestRandomSurv(unittest.TestCase):
    def test_params(self):
        pop_size = 10
        MyNewPathogen = inf.SARS_COV_2(6)   
        sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 3, verbose = 1, enable_surveillance = True, enable_random_testing = True)
        self.assertFalse(sim.pars['enable_random_testing'])     

if __name__ == '__main__':
    unittest.main()