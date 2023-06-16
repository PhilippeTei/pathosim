
import unittest
import infection as inf 
 

class test_diseaseMortalityTests_COVID(unittest.TestCase):
   
    def test_default_death_prob_one(self):
        """
        Infect lots of people with cfr one and short time to die
        duration. Verify that everyone dies, no recoveries.
        """
        pop_size = 200
        pathogen = inf.SARS_COV_2(pop_size)
        n_days = 90 
        sim = inf.Sim(pop_size=pop_size, n_days=n_days, verbose = 0, pathogens = [pathogen]) 

        sim.pathogens[0].rel_symp_prob = 1e6
        sim.pathogens[0].rel_severe_prob = 1e6
        sim.pathogens[0].rel_crit_prob = 1e6
        sim.pathogens[0].rel_death_prob = 1e6
         
        sim.run()
        #print(sim.people['n_infections'])
        assert sim.summary.cum_deaths == pop_size
         

    def test_default_death_prob_zero(self):
        """
        Infect lots of people with probs scaled by 0.2
        Should get 0 deaths, 22 severe cases and 2 critical cases
        """
        pop_size = 10000 
        pathogen = inf.SARS_COV_2(pop_size)
        n_days = 10
        pars = dict(pop_size = pop_size, 
                    n_days = n_days ) 
         
        #run simulation
        sim = inf.Sim(pars=pars, verbose = 0, pathogens = [pathogen]) 
        
        sim.pathogens[0].rel_symp_prob = 0.2
        sim.pathogens[0].rel_severe_prob =0.2
        sim.pathogens[0].rel_crit_prob = 0.2
        sim.pathogens[0].rel_death_prob = 0.2

        sim.run() 

        dead_expected = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        severe_expected = [0., 0., 0., 0., 0., 0., 1., 3., 4., 5., 7.]

        self.assertEqual(True if sim.results['n_severe'] == severe_expected else False, True)
        self.assertEqual(True if sim.results['n_dead'] == dead_expected else False, True) 
         
  
             

if __name__ == '__main__':
    unittest.main()
