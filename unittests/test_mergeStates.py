import unittest 
import infection as inf 
import sciris as sc
import os
import numpy as np 
 
class test_mergeStates(unittest.TestCase):

    def test_mergeStates(self):

        #simple sim to init object
        pop_size = 2
        n_days = 5
        sim = inf.Sim(pop_size=pop_size, pop_infected=0, n_days=n_days, verbose = 0, n_pathogens = 3)
        sim.init_people()

        p1arr1=np.array([True, False]) 
        p2arr1=np.array([False, False]) 
        p3arr1=np.array([False, False]) 

        p1arr2=np.array([False, True]) 
        p2arr2=np.array([False, True]) 
        p3arr2=np.array([False, True]) 

        expectedArr1 = np.array([True, False]) 
        expectedArr2 = np.array([False, True]) 

        sim.people.p_susceptible[0]=p1arr1
        sim.people.p_susceptible[1]=p2arr1
        sim.people.p_susceptible[2]=p3arr1

        
        sim.people.p_recovered[0]=p1arr2
        sim.people.p_recovered[1]=p2arr2
        sim.people.p_recovered[2]=p3arr2

        sim.people.merge_states(3, True, None, False);
         
        self.assertEqual(True if np.array_equal(expectedArr1, sim.people.susceptible)else False, True)
        self.assertEqual(True if np.array_equal(expectedArr2, sim.people.recovered)else False, True)

      
                 
         

if __name__ == '__main__':  
       unittest.main()