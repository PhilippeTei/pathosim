import unittest 
import pathosim as inf 
import sciris as sc
import os
import numpy as np 
 
class test_mergeStates(unittest.TestCase):

    def test_mergeStates(self):

        #simple sim to init object
        pop_size = 2
        n_days = 5

        p1 = inf.Pathogen(1)
        p2 = inf.Pathogen(1)
        p3 = inf.Pathogen(1)

        sim = inf.Sim(pop_size=pop_size, n_days=n_days, verbose = 0, pathogens = [p1,p2,p3])
        sim.initialize()

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

        sim.people.merge_states(3, True, False);
         
        self.assertEqual(True if np.array_equal(expectedArr1, sim.people.susceptible)else False, True) 

      
                 
         

if __name__ == '__main__':  
       unittest.main()